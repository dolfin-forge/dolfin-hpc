// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher 2013. (infinite loop fix)
//
// First added:  2006-05-12
// Last changed: 2013-07-19

#include <dolfin/mesh/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(Mesh& mesh, uint dim) :
    entity_(mesh, dim, 0),
    pos_(0),
    pos_end_(mesh.size(dim)),
    index_(0)
{
  // Compute entities if empty and if number of cells is not zero
  if ((pos_end_ == 0) && (mesh.numCells() != 0))
  {
    pos_end_ = mesh.init(dim);
  }
  // In case we refine the mesh and renumber we lose the mesh entities
  // of the other processes if we do not reinit.

}

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntity& entity, uint dim) :
    entity_(entity.mesh(), dim, 0),
    pos_(0)
{
  // Get connectivity
  MeshConnectivity& c = entity.mesh().topology()(entity.dim(), dim);

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
    pos_end_ = 0;
    index_ = 0;
  }
  else
  {
    pos_end_ = c.size(entity.index());
    index_ = c(entity.index());
  }
}

//-----------------------------------------------------------------------------
MeshEntityIterator::~MeshEntityIterator()
{
}

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntityIterator& entity) :
    entity_(entity.entity_.mesh(), 0, 0),
    pos_(0),
    pos_end_(0),
    index_(0)
{
  error("Illegal use of mesh entity iterator.");
}

//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, const MeshEntityIterator& it)
{
  stream << "[ Mesh entity iterator at position " << it.pos_
      << " stepping from 0 to " << it.pos_end_ - 1 << " ]";
  return stream;
}

//-----------------------------------------------------------------------------

}
