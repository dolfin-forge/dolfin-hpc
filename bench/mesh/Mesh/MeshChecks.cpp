// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#include "MeshChecks.h"

namespace dolfin
{

//-----------------------------------------------------------------------------
bool ghosted_entity_check(MeshEntity& e, bool throw_error)
{
  bool ret = true;
  uint const index = e.index();
  MeshDistributedData const& distdata = e.mesh().distdata();

  // A ghosted entity should have at least one adjacent, the ghost owner
  if (!distdata.is_ghost(e))
  {
    ret = false;
    if (throw_error)
    {
      error("Entity with index %d is not ghosted.", index);
    }
  }

  // A ghosted entity should have at least one adjacent, the ghost owner
  _set<uint> const& adj = distdata.get_shared_adj(e);
  if (adj.count(distdata.get_owner(e)) == 0)
  {
    ret = false;
    if (throw_error)
    {
      error("Ghosted entity with index %d, owner not in adjacent set.", index);
    }
  }

  return ret;
}

//-----------------------------------------------------------------------------
bool shared_entity_check(MeshEntity& e, bool throw_error)
{
  bool ret = true;
  uint const index = e.index();
  MeshDistributedData const& distdata = e.mesh().distdata();

  // A shared entity is... well... shared...
  if (!distdata.is_shared(e))
  {
    ret = false;
    if (throw_error)
    {
      error("Entity with index %d is not shared.", index);
    }
  }

  // A shared entity has a non-empty shared adjacent set
  _set<uint> const& adj = distdata.get_shared_adj(e);
  if (adj.empty())
  {
    ret = false;
    if (throw_error)
    {
      error("Shared entity with index %d has empty adjacent set.", index);
    }
  }

  // A ghosted entity should have at least one adjacent, the ghost owner
  if (distdata.is_ghost(e))
  {
    ghosted_entity_check(e, false);
  }

  return ret;
}

//-----------------------------------------------------------------------------
bool ghosted_entities_check(Mesh& mesh, uint dim, bool throw_error)
{
  message("Check: Ghosted entities of dimension %d.", dim);
  mesh.init(dim);
  // Check consistency of the ghost list
  Array<uint> invalid;
  for (MeshGhostIterator gh(mesh.distdata(), dim); !gh.end(); ++gh)
  {
    MeshEntity e(mesh, dim, gh.index());
    if (!ghosted_entity_check(e, throw_error))
    {
      invalid.push_back(e.index());
    }
  }
  return invalid.empty();
}

//-----------------------------------------------------------------------------
bool shared_entities_check(Mesh& mesh, uint dim, bool throw_error)
{
  message("Check: Shared entities of dimension %d.", dim);
  mesh.init(dim);
  // Test all the mesh entities at the interior boundary
  Array<uint> invalid;
  for (MeshSharedIterator sh(mesh.distdata(), dim); !sh.end(); ++sh)
  {
    MeshEntity e(mesh, dim, sh.index());
    if (!shared_entity_check(e, throw_error))
    {
      invalid.push_back(e.index());
    }
  }
  return invalid.empty();
}

//-----------------------------------------------------------------------------
void append_shared_adj(MeshEntity const& e, std::vector<uint> shared_entities[])
{
  _set<uint> const& adjs = e.mesh().distdata().get_shared_adj(e);
  uint const glob_index = e.mesh().distdata().get_global(e);
  for(_set<uint>::const_iterator s = adjs.begin(); s != adjs.end(); ++s)
  {
    shared_entities[*s].push_back(glob_index);
  }
}

//-----------------------------------------------------------------------------
bool interior_boundary_entities_check(Mesh& mesh, uint dim, bool throw_error)
{
  message("Check: Interior boundary entities of dimension %d.", dim);
  mesh.init(dim);
  BoundaryMesh boundary(mesh, BoundaryMesh::interior);

  // Test all the mesh entities at the interior boundary
  // All the entities should be shared and some are ghosted
  uint const pe_rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();
  std::vector<uint> * shared_entities = new std::vector<uint>[pe_size];

  Array<uint> invalid;
  if (dim == boundary.topology().dim())
  {
    for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
    {
      Facet facet(mesh, boundary.facet_index(*bcell));
      if (!shared_entity_check(facet, throw_error))
      {
        invalid.push_back(facet.index());
      }
      append_shared_adj(facet, shared_entities);
    }
  }
  else
  {
    for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
    {
      uint globIdx = boundary.facet_index(*bcell);
      Facet facet(mesh, globIdx);
      for (MeshEntityIterator e(facet, dim); !e.end(); ++e)
      {
        if (!shared_entity_check(*e, throw_error))
        {
          invalid.push_back(e->index());
        }
        append_shared_adj(*e, shared_entities);
      }
    }
  }

  // Exchange entities and check that the process has the global index

  // Clean
  delete[] shared_entities;

  if (!invalid.empty())
  {
    warning("Interior boundary entities of dim %d (%d) are invalid.", dim,
            invalid.size());
  }
  return invalid.empty();
}

//-----------------------------------------------------------------------------
bool interior_entities_check(Mesh& mesh, uint dim, bool throw_error)
{
  message("Check: Interior entities of dimension %d.", dim);
  error("To be implemented");
  mesh.init(dim);
  BoundaryMesh boundary(mesh, BoundaryMesh::full);

  // Copy indices in a set.

  // Test all the mesh entities at the interior boundary
  Array<uint> invalid;
  return true;
}

//-----------------------------------------------------------------------------

}
