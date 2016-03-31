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
    tdim_(dim),
    gdim_(mesh.geometry().dim()),
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
  uint const * entities = mesh_.topology()(tdim_, entity.tdim_)(index_);
  uint const num_entities = mesh_.topology()(tdim_, entity.tdim_).size(index_);

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
  uint const * entities = mesh_.topology()(tdim_, entity.tdim_)(index_);
  uint const num_entities = mesh_.topology()(tdim_, entity.tdim_).size(index_);

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
uint MeshEntity::global_index() const
{
  return (mesh_.is_distributed() ? mesh_.distdata().get_global(index_, tdim_)
                                 : index_);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_owned() const
{
  return (mesh_.is_distributed() ?
          !mesh_.distdata().is_ghost(index_, tdim_) : true);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_shared() const
{
  return (mesh_.is_distributed() ?
          mesh_.distdata().is_shared(index_, tdim_) : false);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_ghost() const
{
  return (mesh_.is_distributed() ?
          mesh_.distdata().is_ghost(index_, tdim_) : false);
}
//-----------------------------------------------------------------------------
uint MeshEntity::owner() const
{
  return (mesh_.is_distributed() ?
          mesh_.distdata().get_owner(index_, tdim_) : MPI::processNumber());
}
//-----------------------------------------------------------------------------
bool MeshEntity::has_all_vertices_shared() const
{
  if(mesh_.is_distributed())
  {
    if (tdim_ == 0)
    {
      return mesh_.distdata().is_shared(index_, tdim_);
    }
    else
    {
      MeshConnectivity const& c = mesh_.topology()(tdim_, 0);
      for (uint v = 0; v < c.size(index_); ++v)
      {
        if (!mesh_.distdata().is_shared(c(index_)[v], 0))
        {
          return false;
        }
      }
      return true;
    }
  }
  else
  {
    return false;
  }
}
//-----------------------------------------------------------------------------
void MeshEntity::disp() const
{
  section("MeshEntity");
  message("topological dimension : ", tdim_);
  message("geometric dimension   : ", gdim_);
  message("index                 : ", index_);
  begin(  "connectivities        : ");
  for (uint d =0; d < tdim_; ++d)
  {
    cout << d << ": ";
    if(mesh_.topology().is_computed(tdim_, d))
    {
      uint const * entities = mesh_.topology()(tdim_, d)(index_);
      uint const size = mesh_.topology()(tdim_, d).size(index_);
      for (uint i = 0; i < size; ++i)
      {
        cout << "\n\t" << entities[i];
        if(d == 0)
        {
          continue;
        }
        uint const * verts = mesh_.topology()(d, 0)(entities[i]);
        uint const vsize = mesh_.topology()(d, 0).size(entities[i]);
        cout << " ( ";
        for (uint v = 0; v < vsize; ++v)
        {
          cout << verts[v] << ", ";
        }
        cout << ")";
      }
      cout << endl;
    }
    else
    {
      cout << "not computed";
    }
    cout << endl;
  }
  end();
  end();
}
//-----------------------------------------------------------------------------
void MeshEntity::check() const
{
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
