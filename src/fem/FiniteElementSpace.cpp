// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/FiniteElementSpace.h>

#include <dolfin/fem/Form.h>
#include <dolfin/fem/SubSystem.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( Form & form, size_t const i )
  : mesh_( form.dofmaps()[i].mesh() )
  , cell_( mesh_, 0 )
  , finite_element_( new FiniteElement( mesh_.type(), form, i ) )
  , dof_map_( DofMap::acquire( mesh_, form, i ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( Mesh &       mesh,
                                        Form &       form,
                                        size_t const i )
  : mesh_( mesh )
  , cell_( mesh, 0 )
  , finite_element_( new FiniteElement( mesh.type(), form, i ) )
  , dof_map_( DofMap::acquire( mesh, form, i ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( Mesh &                      mesh,
                                        ufc::finite_element const * element,
                                        ufc::dofmap &               dofmap,
                                        bool                        owner )
  : mesh_( mesh )
  , cell_( mesh, 0 )
  , finite_element_( new FiniteElement( *element, owner ) )
  , dof_map_( DofMap::acquire( mesh, dofmap, owner ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( FiniteElementSpace const & space,
                                        size_t const               i )
  : mesh_( space.mesh() )
  , cell_( space.cell() )
  , finite_element_( new FiniteElement( space.element()(), i ) )
  , dof_map_( DofMap::acquire( space.mesh(),
                               *space.dofmap()().create_sub_dofmap( i ),
                               true ) )
{
}

//-----------------------------------------------------------------------------

FiniteElementSpace::FiniteElementSpace( FiniteElementSpace const & space,
                                        SubSystem const &          sub )
  : mesh_( space.mesh() )
  , cell_( space.cell() )
  , finite_element_( new FiniteElement( space.element()(), sub ) )
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
  , finite_element_( &space.element() )
  , dof_map_( DofMap::acquire( other_mesh, *space.dofmap()().create(), true ) )

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
  , finite_element_( &other.element() )
  , dof_map_( DofMap::acquire( other.mesh(), *other.dofmap()().create(), true ) )
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
  std::vector< ufc::finite_element const * > flt_fe =
    finite_element_->flatten();
  std::vector< ufc::dofmap const * > flt_dm = dof_map_.flatten();

  dolfin_assert( flt_fe.size() == flt_dm.size() );
  for ( size_t s = 0; s < flt_fe.size(); ++s )
  {
    flt.push_back( new FiniteElementSpace(
      mesh_, flt_fe[s]->create(), *flt_dm[s]->create(), true ) );
  }
  return flt;
}

//-----------------------------------------------------------------------------

auto FiniteElementSpace::disp() const -> void
{
  section( "FiniteElementSpace" );
  prm( "Finite element", element()().signature() );
  prm( "Dof map", dofmap()().signature() );
  end();
}

//-----------------------------------------------------------------------------

} // end namespace dolfin
