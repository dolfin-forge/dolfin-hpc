// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2009-12-30

#include <dolfin/mesh/VertexNormal.h>

#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/EuclideanBasis.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>


#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(VertexNormal& other) :
    mesh_(other.mesh_),
    subdomain_(NULL),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(none)
{
  *this = other;
}

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(Mesh& mesh, Type weight) :
    mesh_(mesh),
    subdomain_(NULL),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(weight)
{
  uint const nsdim = mesh.topology().dim();

  for (uint d = 0; d < nsdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [nsdim]);
    for (uint i = 0; i < nsdim; ++i)
    {
      basis_.back()[i].init(mesh, 0);
    }
  }
  vertex_type_.init(mesh, 0);

  computeNormal(mesh);
}

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(Mesh& mesh, SubDomain const& subdomain, Type weight) :
    mesh_(mesh),
    subdomain_(&subdomain),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(weight)
{
  uint const nsdim = mesh.topology().dim();

  for (uint d = 0; d < nsdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [nsdim]);
    for (uint i = 0; i < nsdim; ++i)
    {
      basis_.back()[i].init(mesh, 0);
    }
  }
  vertex_type_.init(mesh, 0);

  computeNormal(mesh);
}

//-----------------------------------------------------------------------------
VertexNormal::~VertexNormal()
{
  clear();
}

//-----------------------------------------------------------------------------
void VertexNormal::clear()
{
  while (!basis_.empty())
  {
    delete[] basis_.back();
    basis_.pop_back();
  }
}

//-----------------------------------------------------------------------------
VertexNormal& VertexNormal::operator=(VertexNormal& other)
{
  clear();

  uint const gdim = mesh_.geometry().dim();

  // Initialize data structures
  vertex_type_.init(mesh_, 0);
  for (uint d = 0; d < gdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [gdim]);
    for (uint i = 0; i < gdim; ++i)
    {
      basis_.back()[i].init(mesh_, 0);
    }
  }

  // Copy data
  for (VertexIterator v(mesh_); !v.end(); ++v)
  {
    vertex_type_.set(*v, other.vertex_type_.get(*v));
    for (uint d = 0; d < gdim; ++d)
    {
      for (uint i = 0; i < gdim; ++i)
      {
        basis_[d][i].set(*v, other.basis()[d][i].get(*v));
      }
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
void VertexNormal::getFacetData(VertexNormal::Type type, Mesh& mesh,
                                BoundaryMesh& boundary, Vertex& bvertex,
                                Array<real>& normals, Array<real>& weights)
{
  uint const tdim = mesh.topology().dim();
  uint const gdim = mesh.geometry().dim();
  for (CellIterator bcell(bvertex); !bcell.end(); ++bcell)
  {
    Facet facet(mesh, boundary.facet_index(*bcell));

    if(subdomain_ != NULL && !subdomain_->inside(&bcell->midpoint()[0], true))
    {
      continue;
    }
    dolfin_assert(facet.num_entities(tdim) == 1);
    Cell cell(mesh, facet.entities(tdim)[0]);
    uint local_facet = cell.index(facet);
    Point n = cell.normal(local_facet);
    normals.insert(normals.end(), &n[0], &n[0] + gdim);
    switch (type)
      {
      case VertexNormal::none:
      case VertexNormal::unit:
        weights.push_back(1.0);
        break;
      case VertexNormal::facet:
        weights.push_back(bcell->volume());
        break;
      }
  }
}

//-----------------------------------------------------------------------------
struct VertexData
{
  Array<real> facet_normals;
  Array<real> facet_weights;

  VertexData() :
      facet_normals(),
      facet_weights()
  {
  }
};
typedef _map<uint, VertexData *> VertexDataMap;

//-----------------------------------------------------------------------------
void VertexNormal::computeNormal(Mesh& mesh)
{
  message(1, "VertexNormal: Compute normals");

  uint const tdim = mesh.topology().dim();
  uint const gdim = mesh.geometry().dim();
  // Important: make sure facet to cell connectivities are initialized
  BoundaryMesh& boundary = mesh.exterior_boundary();

  VertexDataMap vdmap;
  int rank = dolfin::MPI::rank();
  int pe_size = dolfin::MPI::size();
  Array<uint> * u_sendbuff = new Array<uint> [pe_size];
  Array<real> * r_sendbuff = new Array<real> [pe_size];

  //--- Collect shared data ---------------------------------------------------
  if (mesh.is_distributed())
  {
#ifdef HAVE_MPI

    DistributedData const& ddv = mesh.distdata()[0];

    // Send buffer for
    // - global index of shared vertices
    // - number of neighbouring boundary cells/global facets
    // Send buffer for :
    // - facets normals associated with shared vertices
    // - weights of facets normals
    if (boundary.num_cells() > 0)
    {
      for (VertexIterator bvertex(boundary); !bvertex.end(); ++bvertex)
      {
        uint const loc_id = boundary.vertex_index(*bvertex);
        //
        if (ddv.is_shared(loc_id))
        {
          if (ddv.is_ghost(loc_id))
          {
            Array<real> normals;
            Array<real> weights;
            getFacetData(type_, mesh, boundary, *bvertex, normals, weights);
            uint const owner = ddv.get_owner(loc_id);
            u_sendbuff[owner].push_back(ddv.get_global(loc_id));
            u_sendbuff[owner].push_back(weights.size());
            r_sendbuff[owner].insert(r_sendbuff[owner].end(), normals.begin(),
                                     normals.end());
            r_sendbuff[owner].insert(r_sendbuff[owner].end(), weights.begin(),
                                     weights.end());
          }
          else
          {
            uint const glb_id = ddv.get_global(loc_id);
            dolfin_assert(glb_id < mesh.global_size(0));
            VertexData * data = new VertexData();
            // Do not fill to avoid copy
            // getFacetData(weighting_, mesh, boundary, *bvertex,
            // data->facet_normals, data->facet_weights);
            vdmap.insert(std::pair<uint, VertexData *>(glb_id, data));
          }
        }
      }
    }

    // Exchange data
    MPI_Status status;
    uint src;
    uint dest;
    uint const u_size = 2;
    int u_recvcount = 0;
    int u_sendcount = 0;
    int r_sendcount = 0;
    int u_maxsendcount = 0;
    int u_maxrecvcount = 0;
    int r_maxsendcount = 0;
    int r_maxrecvcount = 0;
    _set<uint> const& adjs = ddv.get_adj_ranks();
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      u_maxsendcount = std::max(u_maxsendcount, int(u_sendbuff[*it].size()));
      r_maxsendcount = std::max(r_maxsendcount, int(r_sendbuff[*it].size()));
    }
    dolfin_assert(u_maxsendcount <= u_size * mesh.topology().num_ghost(0));
    MPI_Allreduce(&u_maxsendcount, &u_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    dolfin_assert(u_maxrecvcount > 0);
    MPI_Allreduce(&r_maxsendcount, &r_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    dolfin_assert(r_maxrecvcount > 0);

    // For each process
    uint * u_recvbuff = new uint[u_maxrecvcount];
    real * r_recvbuff = new real[r_maxrecvcount];
    for (int j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;

      u_sendcount = u_sendbuff[dest].size();
      MPI_Sendrecv(&u_sendbuff[dest][0], u_sendcount, MPI_UNSIGNED, dest, 1,
                   &u_recvbuff[0], u_maxrecvcount, MPI_UNSIGNED, src, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &u_recvcount);

      r_sendcount = r_sendbuff[dest].size();
      MPI_Sendrecv(&r_sendbuff[dest][0], r_sendcount, MPI_DOUBLE, dest, 1,
                   &r_recvbuff[0], r_maxrecvcount, MPI_DOUBLE, src, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);

      real * rptr = &r_recvbuff[0];
      for (int iiu = 0; iiu < u_recvcount; iiu += u_size)
      {
        uint const glb_id = u_recvbuff[iiu];
        dolfin_assert(glb_id < mesh.global_size(0));
        uint const num_nc = u_recvbuff[iiu + 1];

        VertexDataMap::iterator it = vdmap.find(glb_id);
        if (it != vdmap.end())
        {
          // Add corresponding facet normals and weights
          dolfin_assert(it->second != NULL);
          VertexData * vd = it->second;
          vd->facet_normals.insert(vd->facet_normals.end(), rptr,
                                   rptr + gdim * num_nc);
          rptr += gdim * num_nc;
          vd->facet_weights.insert(vd->facet_weights.end(), rptr,
                                   rptr + num_nc);
          rptr += num_nc;
        }
      }

      // Clear for reuse
      u_sendbuff[dest].clear();
      r_sendbuff[dest].clear();
    }
    delete[] u_recvbuff;
    delete[] r_recvbuff;
#endif
  }

  //--- Compute normals -------------------------------------------------------
  real const cosalpha = std::cos(alpha_max_);
  if (boundary.num_cells() > 0)
  {
    // Initialize cartesian basis
    Point B[EuclideanSpace::MAX_DIMENSION];
    for (uint d = 0; d < EuclideanSpace::MAX_DIMENSION; ++d)
    {
      B[d][d] = 1.0;
    }
    //
    bool weighted = (type_ != VertexNormal::none);
    for (VertexIterator bvertex(boundary); !bvertex.end(); ++bvertex)
    {
      uint const local_index = boundary.vertex_index(*bvertex);
      Vertex v(mesh, local_index);
      if (v.is_ghost())
      {
        continue;
      }

      //--- Get facet normals and weights -----------------------------------
      Array<real> N;
      Array<real> W;
      getFacetData(type_, mesh, boundary, *bvertex, N, W);
      if (v.is_shared())
      {
        VertexDataMap::iterator it = vdmap.find(v.global_index());
        VertexData * vd = it->second;
        N.insert(N.end(), vd->facet_normals.begin(), vd->facet_normals.end());
        W.insert(W.end(), vd->facet_weights.begin(), vd->facet_weights.end());
      }
      dolfin_assert(N.size() == W.size() * gdim);

      //--- Compute basis ---------------------------------------------------
      uint nsurf = EuclideanBasis::compute(gdim, B, N, W, cosalpha, weighted);
      uint vtype = std::min(tdim, nsurf);

      //
      vertex_type_.set(local_index, vtype);
      for (uint e = 0; e < gdim; ++e)
      {
        for (uint d = 0; d < gdim; ++d)
        {
          basis_[e][d].set(local_index, B[e][d]);
        }
      }

      if (v.is_shared())
      {
        _set<uint> const& adjs = mesh.distdata()[0].get_shared_adj(v.index());
        for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end();
             ++it)
        {
          u_sendbuff[*it].push_back(v.global_index());
          u_sendbuff[*it].push_back(vtype);
          for (uint e = 0; e < gdim; ++e)
          {
            real * bptr = &B[e][0];
            r_sendbuff[*it].insert(r_sendbuff[*it].end(), bptr, bptr + gdim);
          }
        }
      }

    }
  }

  if (mesh.is_distributed())
  {
#ifdef HAVE_MPI

    DistributedData const& ddv = mesh.distdata()[0];

    // Exchange data
    MPI_Status status;
    uint src;
    uint dest;
    uint const u_size = 2;
    int u_recvcount = 0;
    int u_sendcount = 0;
    int r_sendcount = 0;
    int u_maxsendcount = 0;
    int u_maxrecvcount = 0;
    int r_maxsendcount = 0;
    int r_maxrecvcount = 0;
    _set<uint> const& adjs = ddv.get_adj_ranks();
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      u_maxsendcount = std::max(u_maxsendcount, int(u_sendbuff[*it].size()));
      r_maxsendcount = std::max(r_maxsendcount, int(r_sendbuff[*it].size()));
    }
    MPI_Allreduce(&u_maxsendcount, &u_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    dolfin_assert(u_maxrecvcount > 0);
    MPI_Allreduce(&r_maxsendcount, &r_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    dolfin_assert(r_maxrecvcount > 0);

    // For each process
    uint * u_recvbuff = new uint[u_maxrecvcount];
    real * r_recvbuff = new real[r_maxrecvcount];
    for (int j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;

      u_sendcount = u_sendbuff[dest].size();
      MPI_Sendrecv(&u_sendbuff[dest][0], u_sendcount, MPI_UNSIGNED, dest, 1,
                   &u_recvbuff[0], u_maxrecvcount, MPI_UNSIGNED, src, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &u_recvcount);

      r_sendcount = r_sendbuff[dest].size();
      MPI_Sendrecv(&r_sendbuff[dest][0], r_sendcount, MPI_DOUBLE, dest, 1,
                   &r_recvbuff[0], r_maxrecvcount, MPI_DOUBLE, src, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);

      //
      uint iir = 0;
      for (int iiu = 0; iiu < u_recvcount; iiu += u_size)
      {
        dolfin_assert(ddv.has_global(u_recvbuff[iiu]));
        uint const local_index = ddv.get_local(u_recvbuff[iiu]);
        dolfin_assert(ddv.is_ghost(local_index));
        vertex_type_.set(local_index, u_recvbuff[iiu + 1]);
        for (uint e = 0; e < gdim; ++e)
        {
          for (uint d = 0; d < gdim; ++d)
          {
            basis_[e][d].set(local_index, r_recvbuff[iir]);
            ++iir;
          }
        }
      }
    }
    delete[] u_recvbuff;
    delete[] r_recvbuff;
#endif
  }

  // Cleanup
  delete[] u_sendbuff;
  delete[] r_sendbuff;
  for (VertexDataMap::iterator it = vdmap.begin(); it != vdmap.end(); ++it)
  {
    delete it->second;
  }
  vdmap.clear();

}

//-----------------------------------------------------------------------------

}

