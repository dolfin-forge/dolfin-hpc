// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(Mesh& mesh, uint dim) :
    entity_(mesh, dim, 0),
    pos_(0),
    end_(mesh.size(dim)),
    index_(NULL)
{
}
//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntity& entity, uint dim) :
    entity_(entity.mesh(), dim, 0),
    pos_(0),
    end_(entity.num_entities(dim)),
    index_(entity.entities(dim))
{
}
//-----------------------------------------------------------------------------
MeshEntityIterator::~MeshEntityIterator()
{
}
//-----------------------------------------------------------------------------
MeshEntityIterator::MeshEntityIterator(MeshEntityIterator& other) :
    entity_(other.entity_.mesh(), 0, 0),
    pos_(0),
    end_(0),
    index_(NULL)
{
  error("Illegal use of mesh entity iterator.");
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
