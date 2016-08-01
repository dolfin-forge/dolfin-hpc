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
    distdata_(mesh.is_distributed() ? &mesh.distdata() : NULL),
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
  MeshConnectivity const& mc = mesh_.topology()(tdim_, entity.tdim_);
  dolfin_assert(mc.size() > 0);
  uint const * entities = mc(index_);
  uint const num_entities = mc.size(index_);

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
  MeshConnectivity const& mc = mesh_.topology()(tdim_, entity.tdim_);
  dolfin_assert(mc.size() > 0);
  uint const * entities = mc(index_);
  uint const num_entities = mc.size(index_);

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
  return mesh_.topology().get_global(*this);
}
//-----------------------------------------------------------------------------
void MeshEntity::global_entities(uint dim, uint * indices) const
{
  // Get list of entities for given topological dimension
  if (distdata_ != NULL)
  {
    MeshConnectivity const& mc = mesh_.topology()(tdim_, dim);
    (*distdata_)[dim].get_global(mc.size(index_), mc(index_), indices);
  }
  else
  {
    MeshConnectivity const& mc = mesh_.topology()(tdim_, dim);
    std::copy(mc(index_), mc(index_) + mc.size(index_), indices);
  }
}
//-----------------------------------------------------------------------------
void MeshEntity::global_entities(uint ** indices) const
{
  // Get list of entities for given topological dimension
  if (distdata_ != NULL)
  {
    for (uint d = 0; d < tdim_; ++d)
    {
      MeshConnectivity const& mc = mesh_.topology()(tdim_, d);
      (*distdata_)[d].get_global(mc.size(index_), mc(index_), indices[d]);
    }
    indices[tdim_][0] = (*distdata_)[tdim_].get_global(index_);
  }
  else
  {
    for (uint d = 0; d < tdim_; ++d)
    {
      MeshConnectivity const& mc = mesh_.topology()(tdim_, d);
      std::copy(mc(index_), mc(index_) + mc.size(index_), indices[d]);
    }
    indices[tdim_][0] = index_;
  }
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_owned() const
{
  return mesh_.topology().is_owned(*this);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_shared() const
{
  return mesh_.topology().is_shared(*this);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_ghost() const
{
  return mesh_.topology().is_ghost(*this);
}
//-----------------------------------------------------------------------------
uint MeshEntity::owner() const
{
  return mesh_.topology().get_owner(*this);
}
//-----------------------------------------------------------------------------
bool MeshEntity::has_all_vertices_shared() const
{
  if(distdata_ != NULL)
  {
    if (tdim_ == 0)
    {
      return (*distdata_)[tdim_].is_shared(index_);
    }
    else
    {
      MeshConnectivity const& c = mesh_.topology()(tdim_, 0);
      dolfin_assert(c.size() > 0);
      for (uint v = 0; v < c.size(index_); ++v)
      {
        if (!(*distdata_)[0].is_shared(c(index_)[v]))
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
bool MeshEntity::on_boundary() const
{
  uint const mdim = mesh_.topology().dim();
  if(distdata_ != NULL)
  {
    if (tdim_ ==  mdim - 1)
    {
      // Facet has one adjacent cell and is not shared, thus is global
      return (this->num_entities(mdim) == 1)
          && !(*distdata_)[tdim_].is_shared(index_);
    }
    else
    {
      MeshConnectivity const& cef = mesh_.topology()(tdim_, mdim - 1);
      MeshConnectivity const& cfc = mesh_.topology()(mdim - 1, mdim);
      for (uint f = 0; f < cef.size(index_); ++f)
      {
        uint const fidx = cef(index_)[f];
        if ((cfc.size(fidx) == 1) && !(*distdata_)[mdim - 1].is_shared(fidx))
        {
          // Facet has one adjacent cell and is not shared, thus is global
          return true;
        }
      }
    }
  }
  else
  {
    if (tdim_ == mdim - 1)
    {
      // Facet has one adjacent cell only, thus is global
      return (this->num_entities(mdim) == 1);
    }
    else
    {
      MeshConnectivity const& cef = mesh_.topology()(tdim_, mdim - 1);
      MeshConnectivity const& cfc = mesh_.topology()(mdim - 1, mdim);
      for (uint f = 0; f < cef.size(index_); ++f)
      {
        uint const fidx = cef(index_)[f];
        if ((cfc.size(fidx) == 1))
        {
          // Facet has one adjacent cell only
          return true;
        }
      }
    }
  }
  return false;
}
//-----------------------------------------------------------------------------
void MeshEntity::disp() const
{
  section("MeshEntity");
  //---
  message("topological dimension : %u", tdim_);
  message("geometric dimension   : %u", gdim_);
  message("index                 : %u", index_);
  begin(  "connectivities        : %u");
  MeshTopology const& topology = mesh_.topology();
  for (uint d =0; d < tdim_; ++d)
  {
    cout << d << ": ";
    if(topology.is_computed(tdim_, d))
    {
      uint const * entities = topology(tdim_, d)(index_);
      uint const size = topology(tdim_, d).size(index_);
      for (uint i = 0; i < size; ++i)
      {
        cout << "\n\t" << entities[i];
        if(d == 0)
        {
          continue;
        }
        uint const * verts = topology(d, 0)(entities[i]);
        uint const vsize = topology(d, 0).size(entities[i]);
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
  //---
  end();
}
//-----------------------------------------------------------------------------
void MeshEntity::check() const
{
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
