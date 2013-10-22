// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-12
// Last changed: 2007-04-12

#include <dolfin/mesh/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(Mesh& mesh, uint dim) :
    entity(mesh, dim, 0),
    _pos(0),
    pos_end(mesh.size(dim)),
    index(0)
{
  // FIXME: NOT GOOD
  // Compute entities if empty and if number of cells is not zero
  if ((pos_end == 0) && (mesh.numCells() != 0))
  {
    pos_end = mesh.init(dim);
  }
  // In case we refine the mesh and renumber we lose the mesh entities
  // of the other processes if we do not reinit.

}

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntity& entity, uint dim) :
    entity(entity.mesh(), dim, 0),
    _pos(0)
{
  // Get connectivity
  MeshConnectivity& c = entity.mesh().topology()(entity.dim(), dim);

  // FIXME: NOT GOOD
  // Compute connectivity if empty
  if (c.size() == 0)
  {
    entity.mesh().init(entity.dim(), dim);
  }
  // In case we refine the mesh and renumber we lose the mesh entities
  // if we do not reinit.

  // Get size and index map
  if (c.size() == 0)
  {
    pos_end = 0;
    index = 0;
  }
  else
  {
    pos_end = c.size(entity.index());
    index = c(entity.index());
  }
}

//-----------------------------------------------------------------------------
MeshEntityIterator::~MeshEntityIterator()
{
}

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntityIterator& entity) :
    entity(entity.entity.mesh(), 0, 0),
    _pos(0),
    pos_end(0),
    index(0)
{
  error("Illegal use of mesh entity iterator.");
}

//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, const MeshEntityIterator& it)
{
  stream << "[ Mesh entity iterator at position " << it._pos
      << " stepping from 0 to " << it.pos_end - 1 << " ]";
  return stream;
}

//-----------------------------------------------------------------------------

}
