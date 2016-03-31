// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#include "MeshChecks.h"

namespace dolfin
{

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
  uint const tdim = mesh.topology().dim();
  BoundaryMesh boundary(mesh, BoundaryMesh::interior);
  Array<uint> invalid_shared;
  Array<uint> invalid_neighb;
  if(boundary.numCells() > 0)
  {
    // A bug causes segmentation fault if the boundary is empty
    boundary.init(dim);
    uint const bdim = boundary.topology().dim();
    if (dim > bdim)
    {
      error("interior_boundary_entities_check only works for facets.");
    }

    uint const num_shared = mesh.topology().num_shared(dim);
    uint const num_intbnd = boundary.topology().num_local(dim);
    if(num_shared != num_intbnd)
    {
      error("Inconsistent number of entities: (shared) %d  != %d (boundary)",
            num_shared, num_intbnd);
    }

    // Test all the mesh entities at the interior boundary
    // All the entities should be shared and some are ghosted
    MeshDistributedData& distdata = mesh.distdata();
    if (dim == boundary.topology().dim())
    {
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        if (!distdata.check_shared(facet.index(), facet.dim(), throw_error))
        {
          invalid_shared.push_back(facet.index());
        }
        if (facet.num_entities(tdim) != 1)
        {
          invalid_neighb.push_back(facet.index());
        }
      }
    }
    else
    {
      mesh.init(boundary.topology().dim(), dim);
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        for (MeshEntityIterator e(facet, dim); !e.end(); ++e)
        {
          if (!distdata.check_shared(e->index(), e->dim(), throw_error))
          {
            invalid_shared.push_back(e->index());
          };
        }
      }
    }

    if (!invalid_shared.empty())
    {
      warning("Interior boundary entities of dim %d: %d invalid shared data.",
              dim, invalid_shared.size());
    }
    if (!invalid_neighb.empty())
    {
      warning("Interior boundary entities of dim %d: %d invalid connectivity.",
              dim, invalid_shared.size());
    }
  }
  return invalid_shared.empty() && invalid_neighb.empty();
}

//-----------------------------------------------------------------------------
bool exterior_boundary_entities_check(Mesh& mesh, uint dim, bool throw_error)
{
  message("Check: Exterior boundary entities of dimension %d.", dim);
  mesh.init(dim);
  BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
  Array<uint> invalid;
  if(boundary.numCells() > 0)
  {
    // A bug causes segmentation fault if the boundary is empty
    boundary.init(dim);

    // Test all the mesh entities at the interior boundary
    // All the entities should be shared and some are ghosted
    if (dim == boundary.topology().dim())
    {
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        if(facet.is_shared())
        {
          invalid.push_back(facet.index());
        }
      }
    }
    else
    {
      mesh.init(boundary.topology().dim(), dim);
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
      }
    }

    if (!invalid.empty())
    {
      warning("Exterior boundary entities of dim %d (%d) are invalid.", dim,
              invalid.size());
    }
  }
  return invalid.empty();
}

//-----------------------------------------------------------------------------

}
