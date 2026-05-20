// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/FiniteElementSpace.h>

#include <dolfin/fem/Form.h>
#include <dolfin/fem/SubSystem.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------


FiniteElementSpace::FiniteElementSpace( FiniteElement const & element,
                                        DofMap &              dof )
  : mesh_( dof.mesh() )
  , cell_( dof.mesh(), 0 )
  , finite_element_( std::make_shared<FiniteElement>( element ) )
  , dof_map_( DofMap::acquire( mesh_, dof.ufc(), false ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( Mesh &                      mesh,
                                        ufc::finite_element const * element,
                                        ufc::dofmap const &         dofmap )
  : mesh_( mesh )
  , cell_( mesh, 0 )
  , finite_element_( std::make_shared<FiniteElement>( *element ) )
  , dof_map_( DofMap::acquire( mesh, dofmap, false ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( FiniteElementSpace const & space,
                                        size_t const               i )
  : mesh_( space.mesh() )
  , cell_( space.cell() )
  , finite_element_( std::make_shared<FiniteElement>( space.element().ufc(), i ) )
  , dof_map_( DofMap::acquire( space.mesh(),
                               *space.dofmap().ufc().create_sub_dofmap( i ),
                               true ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( FiniteElementSpace const & space,
                                        SubSystem const &          sub )
  : mesh_( space.mesh() )
  , cell_( space.cell() )
  , finite_element_( std::make_shared<FiniteElement>( space.element().ufc(), sub ) )
  , dof_map_( DofMap::acquire( space.mesh(),
                               *space.dofmap().create_sub_dofmap( sub ),
                               true ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( Mesh &                     other_mesh,
                                        FiniteElementSpace const & space )
  : mesh_( other_mesh )
  , cell_( space.cell() )
  , finite_element_( space.element_ptr() )
  , dof_map_( DofMap::acquire( other_mesh, space.dofmap().ufc(), false ) )
{
  if ( other_mesh.type().cellType() != space.cell().type() )
  {
    error( "Provided mesh is incompatible with given space" );
  }
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( FiniteElementSpace const & other )
  : mesh_( other.mesh() )
  , cell_( other.cell() )
  , finite_element_( other.element_ptr() )
  , dof_map_( DofMap::acquire( other.mesh(), other.dofmap().ufc(), false ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::~FiniteElementSpace()
{
  DofMap::release( dof_map_ );
}

//-----------------------------------------------------------------------------

auto FiniteElementSpace::flatten() const -> std::vector< FiniteElementSpace * >
{
  std::vector< FiniteElementSpace * >        flt;
  std::vector< ufc::finite_element const * > flt_fe = finite_element_->flatten();
  std::vector< ufc::dofmap const * >         flt_dm = dof_map_.flatten();

  dolfin_assert( flt_fe.size() == flt_dm.size() );
  for ( size_t s = 0; s < flt_fe.size(); ++s )
  {
    flt.push_back( new FiniteElementSpace( mesh_, flt_fe[s], *flt_dm[s] ) );
  }

  return flt;
}

//-----------------------------------------------------------------------------

auto FiniteElementSpace::disp() const -> void
{
  section( "FiniteElementSpace" );
  prm( "Finite element", element().signature );
  prm( "Dof map", dofmap().signature );
  end();
}

//-----------------------------------------------------------------------------

} // end namespace dolfin
