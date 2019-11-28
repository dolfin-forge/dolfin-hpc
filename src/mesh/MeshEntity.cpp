// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEntity::MeshEntity(Mesh& mesh, uint dim, uint index) :
    mesh_(mesh),
    topology_(mesh_.topology()),
    tdim_(dim),
    gdim_(mesh.geometry_dimension()),
    distdata_(topology_.distdata_),
    index_(index)
{
}
//-----------------------------------------------------------------------------
MeshEntity::~MeshEntity()
{
}
//-----------------------------------------------------------------------------
void MeshEntity::get_entities(uint dim, uint * indices) const
{
  uint *b, *e;
  topology_(tdim_, dim)(index_, b, e);
  std::copy(b, e, indices);
}
//-----------------------------------------------------------------------------
void MeshEntity::get_entities(uint ** indices) const
{
  uint *b, *e;
  for (uint dim = 0; dim < tdim_; ++dim)
  {
    topology_(tdim_, dim)(index_, b, e);
    std::copy(b, e, indices[dim]);
  }
  indices[tdim_][0] = index_;
}
//-----------------------------------------------------------------------------
bool MeshEntity::incident(MeshEntity const& entity) const
{
  // Must be in the same mesh to be incident
  if (&topology_ != &entity.topology_) return false;

  return topology_(tdim_, entity.tdim_).incident(index_, entity.index_);
}
//-----------------------------------------------------------------------------
int MeshEntity::index(MeshEntity const& entity) const
{
  // Must be in the same mesh to be incident
  if (&topology_ != &entity.topology_)
  {
    error("Unable to compute index of an entity defined on a different mesh.");
  }

  return topology_(tdim_, entity.tdim_).index(index_, entity.index_);
}
//-----------------------------------------------------------------------------
uint MeshEntity::global_index() const
{
  return (topology_.distributed() ? distdata_[tdim_].get_global(index_) : index_);
}
//-----------------------------------------------------------------------------
void MeshEntity::get_global_entities(uint dim, uint * indices) const
{
  // Get list of entities for given topological dimension
  if ( topology_.distributed() )
  {
    Connectivity const& mc = topology_(tdim_, dim);
    distdata_[dim].get_global(mc.degree(index_), mc(index_), indices);
  }
  else
  {
    get_entities(dim, indices);
  }
}
//-----------------------------------------------------------------------------
void MeshEntity::get_global_entities(uint ** indices) const
{
  // Get list of entities for given topological dimension
  if ( topology_.distributed() )
  {
    for (uint d = 0; d < tdim_; ++d)
    {
      Connectivity const& mc = topology_(tdim_, d);
      distdata_[d].get_global(mc.degree(index_), mc(index_), indices[d]);
    }
    indices[tdim_][0] = distdata_[tdim_].get_global(index_);
  }
  else
  {
    get_entities(indices);
  }
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_owned() const
{
  return (topology_.distributed() ? distdata_[tdim_].is_owned(index_) : true);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_shared() const
{
  return (topology_.distributed() ? distdata_[tdim_].is_shared(index_) : false);
}
//-----------------------------------------------------------------------------
bool MeshEntity::is_ghost() const
{
  return (topology_.distributed() ? distdata_[tdim_].is_ghost(index_) : false);
}
//-----------------------------------------------------------------------------
uint MeshEntity::owner() const
{
  return (topology_.distributed() ? distdata_[tdim_].get_owner(index_) : MPI::rank());
}
//-----------------------------------------------------------------------------
_set<uint> const * MeshEntity::adjacents() const
{
  return (topology_.distributed() ? distdata_[tdim_].ptr_shared_adj(index_) : NULL);
}
//-----------------------------------------------------------------------------
bool MeshEntity::has_all_vertices_shared() const
{
  if( topology_.distributed() )
  {
    if (tdim_ == 0)
    {
      return distdata_[tdim_].is_shared(index_);
    }
    else
    {
      Connectivity const& c = topology_(tdim_, 0);
      dolfin_assert(c.order() > 0);
      for (uint v = 0; v < c.degree(index_); ++v)
      {
        if ( not distdata_[0].is_shared(c(index_)[v]) )
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
  uint const mdim = topology_.dim();
  uint const fdim = topology_.type(index_).facet_dim();
  if( topology_.distributed() )
  {
    if (tdim_ ==  fdim)
    {
      // Facet has one adjacent cell and is not shared, thus is global
      return (this->num_entities(mdim) == 1)
          && (not distdata_[tdim_].is_shared(index_));
    }
    else
    {
      Connectivity const& cef = topology_(tdim_, fdim);
      Connectivity const& cfc = topology_(fdim , mdim);
      for (uint f = 0; f < cef.degree(index_); ++f)
      {
        uint const fidx = cef(index_)[f];
        if ( (cfc.degree(fidx) == 1) && (not distdata_[fdim].is_shared(fidx)) )
        {
          // Facet has one adjacent cell and is not shared, thus is global
          return true;
        }
      }
    }
  }
  else
  {
    if (tdim_ == fdim)
    {
      // Facet has one adjacent cell only, thus is global
      return (this->num_entities(mdim) == 1);
    }
    else
    {

      Connectivity const& cef = topology_(tdim_, fdim);
      Connectivity const& cfc = topology_(fdim , mdim);
      for (uint f = 0; f < cef.degree(index_); ++f)
      {
        uint const fidx = cef(index_)[f];
        if ((cfc.degree(fidx) == 1))
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
  prm("topological dimension" , tdim_);
  prm("geometric dimension"   , gdim_);
  prm("index"                 , index_);
  MeshTopology const& topology = topology_;
  for (uint d =0; d < tdim_; ++d)
  {
    cout << d << ": ";
    if(topology.connectivity(tdim_, d))
    {
      uint const * entities = topology(tdim_, d)(index_);
      uint const size = topology(tdim_, d).degree(index_);
      for (uint i = 0; i < size; ++i)
      {
        cout << "\n\t" << entities[i];
        if(d == 0)
        {
          continue;
        }
        uint const * verts = topology(d, 0)(entities[i]);
        uint const vsize = topology(d, 0).degree(entities[i]);
        cout << " ( ";
        for (uint v = 0; v < vsize; ++v)
        {
          cout << verts[v] << ", ";
        }
        cout << ")";
      }
      cout << "\n";
    }
    else
    {
      cout << "not computed";
    }
    cout << "\n";
  }
  end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
