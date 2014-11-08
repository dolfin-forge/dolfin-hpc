// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-11
// Last changed: 2006-10-20

#include <dolfin/mesh/MeshEntity.h>

#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEntity::MeshEntity(Mesh& mesh, uint dim, uint index) :
    mesh_(mesh),
    dim_(dim),
    index_(index)
{
}
//-----------------------------------------------------------------------------
MeshEntity::~MeshEntity()
{
}
//-----------------------------------------------------------------------------
bool MeshEntity::incident(MeshEntity const& entity) const
{
  // Must be in the same mesh to be incident
  if (&mesh_ != &entity.mesh_) return false;

  // Get list of entities for given topological dimension
  uint const * entities = mesh_.topology()(dim_, entity.dim_)(index_);
  uint const num_entities = mesh_.topology()(dim_, entity.dim_).size(index_);

  // Check if any entity matches
  for (uint i = 0; i < num_entities; ++i)
  {
    if (entities[i] == entity.index_) return true;
  }

  // Entity was not found
  return false;
}
//-----------------------------------------------------------------------------
uint MeshEntity::index(MeshEntity const& entity) const
{
  // Must be in the same mesh to be incident
  if (&mesh_ != &entity.mesh_)
  {
    error("Unable to compute index of given entity defined on a different "
          "mesh.");
  }

  // Get list of entities for given topological dimension
  uint const * entities = mesh_.topology()(dim_, entity.dim_)(index_);
  uint const num_entities = mesh_.topology()(dim_, entity.dim_).size(index_);

  // Check if any entity matches
  for (uint i = 0; i < num_entities; ++i)
  {
    if (entities[i] == entity.index_) return i;
  }

  // Entity was not found
  error("Unable to compute index of given entity (not found).");

  return 0;
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_shared() const
{
  return mesh_.distdata().is_shared(index_, dim_);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_ghost() const
{
  return mesh_.distdata().is_ghost(index_, dim_);
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, MeshEntity const& entity)
{
  stream << "[ Mesh entity " << entity.index() << " of topological dimension "
         << entity.dim() << " ]";
  return stream;
}
//-----------------------------------------------------------------------------

}
