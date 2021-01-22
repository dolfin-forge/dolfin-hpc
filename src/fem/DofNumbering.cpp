// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/mesh/Mesh.h>

// For factory function
#include <dolfin/fem/numbering/CG1s.h>
#include <dolfin/fem/numbering/CG1v.h>
#include <dolfin/fem/numbering/DG0s.h>
#include <dolfin/fem/numbering/DG0v.h>
#include <dolfin/fem/numbering/Parallel0.h>
#include <dolfin/fem/numbering/Parallel1.h>
#include <dolfin/fem/numbering/RealSpace.h>
#include <dolfin/fem/numbering/Serial.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

auto DofNumbering::create(Mesh& mesh, ufc::dofmap& ufc_dofmap) -> DofNumbering *
{
  DofNumbering * ret = nullptr;
  size_t const tdim = mesh.topology_dimension();
  size_t const num_verts = mesh.topology().global_size(0);
  size_t const num_cells = mesh.topology().global_size(tdim);

  std::vector<ufc::dofmap const*> flattened( 1, nullptr );
  flattened.resize( 0 );

  // FIXME num_entities should maybe be stored somewhere else
  std::vector< size_t > num_entities( tdim+1, 0 );
  for ( size_t d = 0; d <= tdim; ++d )
  {
    if ( mesh.topology().connectivity( d ) )
    {
      num_entities[d] = mesh.topology().global_size( d );
    }
  }

  ///
  DofMap::flatten(&ufc_dofmap, flattened);
  size_t const value_size = flattened.size();
  bool const vector = DofMap::can_vectorize(flattened);
  // destruct( flattened );
  for( size_t i = 0; i < flattened.size(); ++i )
    if( flattened[i] != nullptr )
      delete flattened[i];

  // UFC dofmap should be initialized to compute the global dimension
  DofNumbering::init(mesh, ufc_dofmap);

  // Real
  if (ufc_dofmap.global_dimension(num_entities) == ufc_dofmap.num_element_dofs())
  {
    message(1, "DofNumbering: create RealSpaceNumbering for %s",
            ufc_dofmap.signature());
    ret = new RealSpaceNumbering(mesh, ufc_dofmap);
  }
  // CG1s
  else if (ufc_dofmap.global_dimension(num_entities) == num_verts)
  {
    message(1, "DofNumbering: create CG1sNumbering for %s",
            ufc_dofmap.signature());
    ret = new CG1sNumbering(mesh, ufc_dofmap);
  }
  // DG0s
  else if (ufc_dofmap.global_dimension(num_entities) == num_cells)
  {
    message(1, "DofNumbering: create DG0sNumbering for %s",
            ufc_dofmap.signature());
    ret = new DG0sNumbering(mesh, ufc_dofmap);
  }
  // DG0v
  else if (vector && ufc_dofmap.global_dimension(num_entities) == value_size * num_cells)
  {
    message(1, "DofNumbering: create DG0vNumbering for %s",
            ufc_dofmap.signature());
    ret = new DG0vNumbering(mesh, ufc_dofmap);
  }
  else if (mesh.is_distributed())
  {
    // Generic parallel for scalar
    if (value_size == 1)
    {
      message(1, "DofNumbering: create Parallel0Numbering for %s",
              ufc_dofmap.signature());
      ret = new Parallel0Numbering(mesh, ufc_dofmap);
    }
    // CG1v
    else if (vector && ufc_dofmap.global_dimension(num_entities) == value_size * num_verts)
    {
      message(1, "DofNumbering: create CG1vNumbering for %s",
              ufc_dofmap.signature());
      ret = new CG1vNumbering(mesh, ufc_dofmap);
    }
    // Generic parallel for value rank > 0
    else
    {
      message(1, "DofNumbering: create Parallel1Numbering for %s",
              ufc_dofmap.signature());
      ret = new Parallel1Numbering(mesh, ufc_dofmap);
    }
  }
  else
  {
    message(1, "DofNumbering: create SerialNumbering for %s",
            ufc_dofmap.signature());
    ret = new SerialNumbering(mesh, ufc_dofmap);
  }
  return ret;
}

//-----------------------------------------------------------------------------

DofNumbering::DofNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
    mesh(mesh),
    ufc_dofmap(ufc_dofmap),
    array_size(0),
    array(nullptr),
    offset_(0),
    size_(0)
{
}

//-----------------------------------------------------------------------------

DofNumbering::~DofNumbering()
{
  clear();
}

//-----------------------------------------------------------------------------

DofNumbering::DofNumbering(DofNumbering const& other) :
    mesh(other.mesh),
    ufc_dofmap(other.ufc_dofmap),
    array_size(other.array_size),
    array(nullptr),
    offset_(other.offset_),
    size_(other.size_)
{
  if (other.array != nullptr)
  {
    array = new size_t[array_size];
    std::copy(other.array, other.array + other.array_size, array);
  }
}

//-----------------------------------------------------------------------------

void DofNumbering::pretabulate(size_t *& array, size_t& array_size) const
{
  if (this->size_ == 0)
  {
    error("DofNumbering::cache : numbering range is zero");
  }
  array = new size_t[mesh.num_cells() * ufc_dofmap.num_element_dofs()];
  array_size = 0;
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);
  for (; !cell.end(); ++cell, array_size += ufc_dofmap.num_element_dofs())
  {
    ufc_cell.update(*cell);
    tabulate_dofs(&array[array_size], ufc_cell, *cell);
  }
}

//-----------------------------------------------------------------------------

void DofNumbering::disp() const
{
  section("DofNumbering");
  message("Dofmap : %s", ufc_dofmap.signature());
  message("Type   : %s", this->description().c_str());
  end();
}

//-----------------------------------------------------------------------------

void DofNumbering::init(Mesh& mesh, ufc::dofmap& ufc_dofmap)
{
  if (ufc_dofmap.topological_dimension() != mesh.topology_dimension())
  {
    error("DofNumbering::init : invalid topological dimension.");
  }
  // if (ufc_dofmap.geometric_dimension() != mesh.geometry_dimension())
  // {
  //   error("DofNumbering::init : invalid geometrical dimension.");
  // }
  if (ufc_dofmap.topological_dimension() != mesh.geometry_dimension())
  {
    error("DofNumbering::init : topological dimension != geometric dimension.\n"
          "Unfortunately FFC cannot generate that (yeah that sux poneys).");
  }

  // Initialize mesh entities used by dof map
  for (size_t d = 0; d <= mesh.topology_dimension(); ++d)
  {
    if (ufc_dofmap.needs_mesh_entities(d))
    {
      mesh.init(d);
    }
  }

  // FIXME num_entities should maybe be stored somewhere else
  std::vector< size_t > num_entities( ufc_dofmap.topological_dimension()+1, 0 );
  for ( size_t d = 0; d <= ufc_dofmap.topological_dimension(); ++d )
  {
    if ( mesh.topology().connectivity( d ) )
    {
      num_entities[d] = mesh.topology().global_size( d );
    }
  }

  message(1, "DofNumbering: initialized UFC dofmap with global dimension %u",
          ufc_dofmap.global_dimension(num_entities));
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
