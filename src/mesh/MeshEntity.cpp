// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshEntity.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>

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
        if ( not distdata_[0].is_shared(c[index_][v]) )
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
  uint const fdim = topology_.type().facet_dim();
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
        uint const fidx = cef[index_][f];
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
        uint const fidx = cef[index_][f];
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
      Array<uint> const & entities = topology(tdim_, d)[index_];
      uint const size = topology(tdim_, d).degree(index_);
      for (uint i = 0; i < size; ++i)
      {
        cout << "\n\t" << entities[i];
        if(d == 0)
        {
          continue;
        }
        Array<uint> const & verts = topology(d, 0)[entities[i]];
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
