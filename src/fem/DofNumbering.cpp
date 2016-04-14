// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Martin Alnes, 2008
// Modified by Niclas Jansson, 2009
// Modified by Aurélien Larcher, 2013  (Pk bug, extension and partial rewrite)
//
// First added:  2007-03-01
// Last changed: 2009-11-01

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/UFCMesh.h>
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
DofNumbering * DofNumbering::create(Mesh& mesh, ufc::dofmap& ufc_dofmap)
{
  DofNumbering * ret = NULL;
  uint const tdim = mesh.topology().dim();
  uint const num_verts = mesh.topology().global_size(0);
  uint const num_cells = mesh.topology().global_size(tdim);

  ///
  Array<ufc::dofmap const*> flattened;
  DofMap::flatten(&ufc_dofmap, flattened);
  uint const value_size = flattened.size();
  bool const vector = DofMap::can_vectorize(flattened);
  flattened.free();

  // UFC dofmap should be initialized to compute the global dimension
  DofNumbering::init(mesh, ufc_dofmap);

  // Real
  if (ufc_dofmap.global_dimension() == ufc_dofmap.local_dimension())
  {
    message(1, "DofNumbering : create RealSpaceNumbering for %s",
            ufc_dofmap.signature());
    ret = new RealSpaceNumbering(mesh, ufc_dofmap);
  }
  // CG1s
  else if (ufc_dofmap.global_dimension() == num_verts)
  {
    message(1, "DofNumbering : create CG1sNumbering for %s",
            ufc_dofmap.signature());
    ret = new CG1sNumbering(mesh, ufc_dofmap);
  }
  // DG0s
  else if (ufc_dofmap.global_dimension() == num_cells)
  {
    message(1, "DofNumbering : create DG0sNumbering for %s",
            ufc_dofmap.signature());
    ret = new DG0sNumbering(mesh, ufc_dofmap);
  }
  // CG1v
  else if (vector && ufc_dofmap.global_dimension() == value_size * num_verts)
  {
    message(1, "DofNumbering : create CG1vNumbering for %s",
            ufc_dofmap.signature());
    ret = new CG1vNumbering(mesh, ufc_dofmap);
  }
  // DG0v
  else if (vector && ufc_dofmap.global_dimension() == value_size * num_cells)
  {
    message(1, "DofNumbering : create DG0vNumbering for %s",
            ufc_dofmap.signature());
    ret = new DG0vNumbering(mesh, ufc_dofmap);
  }
  else if (mesh.is_distributed())
  {
    // Generic parallel for scalar
    if (value_size == 1)
    {
      message(1, "DofNumbering : create Parallel0Numbering for %s",
              ufc_dofmap.signature());
      ret = new Parallel0Numbering(mesh, ufc_dofmap);
    }
    // Generic parallel for value rank > 0
    else
    {
      message(1, "DofNumbering : create Parallel1Numbering for %s",
              ufc_dofmap.signature());
      ret = new Parallel1Numbering(mesh, ufc_dofmap);
    }
  }
  else
  {
    message(1, "DofNumbering : create SerialNumbering for %s",
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
    array(NULL),
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
    array(NULL),
    offset_(other.offset_),
    size_(other.size_)
{
  if (other.array != NULL)
  {
    array = new uint[array_size];
    std::copy(other.array, other.array + other.array_size, array);
  }
}
//-----------------------------------------------------------------------------
DofNumbering& DofNumbering::operator=(DofNumbering const& other)
{
  return *this;
}
//-----------------------------------------------------------------------------
uint DofNumbering::offset() const
{
  return offset_;
}
//-----------------------------------------------------------------------------
uint DofNumbering::size() const
{
  return size_;
}
//-----------------------------------------------------------------------------
uint const * DofNumbering::block() const
{
  if (array == NULL)
  {
    pretabulate(array, array_size);
  }
  return array;
}
//-----------------------------------------------------------------------------
uint DofNumbering::block_size() const
{
  if (array == NULL)
  {
    pretabulate(array, array_size);
  }
  return array_size;
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
void DofNumbering::init()
{
  DofNumbering::clear();
  DofNumbering::init(mesh, ufc_dofmap);
}
//-----------------------------------------------------------------------------
void DofNumbering::clear()
{
  offset_ = 0;
  size_ = 0;
  array_size = 0;
  delete [] array;
  array = NULL;
}
//-----------------------------------------------------------------------------
void DofNumbering::set_range(uint offset, uint size)
{
  offset_ = offset;
  size_ = size;
}
//-----------------------------------------------------------------------------
void DofNumbering::pretabulate(uint *& array, uint& array_size) const
{
  if (this->size_ == 0)
  {
    error("DofNumbering::cache : numbering range is zero");
  }
  array = new uint[mesh.num_cells() * ufc_dofmap.local_dimension()];
  array_size = 0;
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);
  for (; !cell.end(); ++cell, array_size += ufc_dofmap.local_dimension())
  {
    ufc_cell.update(*cell);
    tabulate_dofs(&array[array_size], ufc_cell, *cell);
  }
}
//-----------------------------------------------------------------------------
void DofNumbering::init(Mesh& mesh, ufc::dofmap& ufc_dofmap)
{
  if (ufc_dofmap.topological_dimension() != mesh.topology().dim())
  {
    error("DofNumbering::init : invalid topological dimension.");
  }
  if (ufc_dofmap.geometric_dimension() != mesh.geometry().dim())
  {
    error("DofNumbering::init : invalid topological dimension.");
  }
  if (ufc_dofmap.topological_dimension() != mesh.geometry().dim())
  {
    error("DofNumbering::init : topological dimension != geometric dimension.\n"
          "Unfortunately FFC cannot generate that (yeah that sux poneys).");
  }

  // Initialize mesh entities used by dof map
  for (uint d = 0; d <= mesh.topology().dim(); ++d)
  {
    if (ufc_dofmap.needs_mesh_entities(d))
    {
      mesh.init(d);
    }
  }

  // Initialize UFC mesh data after entities are created since global sizes may
  // otherwise be incorrect.
  UFCMesh ufc_mesh(mesh);

  // Initialize UFC dof map
  if (ufc_dofmap.init_mesh(ufc_mesh))
  {
    CellIterator cell(mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell);
      ufc_dofmap.init_cell(ufc_mesh, ufc_cell);
    }
    ufc_dofmap.init_cell_finalize();
  }

  message("DofNumbering : initialized UFC dofmap with global dimension %u",
          ufc_dofmap.global_dimension());
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
