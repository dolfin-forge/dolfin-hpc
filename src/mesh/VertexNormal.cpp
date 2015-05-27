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
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>

#define DEBUG 1

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
    dolfin_assert(facet.numEntities(tdim) == 1);
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
  mesh.renumber();

  uint const tdim = mesh.topology().dim();
  uint const gdim = mesh.geometry().dim();
  // Important: make sure facet to cell connectivities are initialized
  BoundaryMesh& boundary = mesh.exterior_boundary();
  mesh.init(boundary.topology().dim(), tdim);

  VertexDataMap vdata;
  int rank = dolfin::MPI::processNumber();
  int pe_size = dolfin::MPI::numProcesses();
  Array<uint> * u_sendbuff = new Array<uint> [pe_size];
  Array<real> * r_sendbuff = new Array<real> [pe_size];

  //--- Collect shared data ---------------------------------------------------
  if (mesh.is_distributed())
  {
#ifdef HAVE_MPI

    // Send buffer for
    // - global index of shared vertices
    // - number of neighbouring boundary cells/global facets
    // Send buffer for :
    // - facets normals associated with shared vertices
    // - weights of facets normals
    if (boundary.numCells() > 0)
    {
      boundary.init(boundary.topology().dim(), 0);
      for (VertexIterator bvertex(boundary); !bvertex.end(); ++bvertex)
      {
        uint const loc_id = boundary.vertex_index(*bvertex);
        //
        if (mesh.distdata().is_shared(loc_id, 0))
        {
          if (mesh.distdata().is_ghost(loc_id, 0))
          {
            Array<real> normals;
            Array<real> weights;
            getFacetData(type_, mesh, boundary, *bvertex, normals, weights);
            uint const glob_id = mesh.distdata().get_vertex_global(loc_id);
            uint const owner = mesh.distdata().get_owner(loc_id, 0);
            u_sendbuff[owner].push_back(glob_id);
            u_sendbuff[owner].push_back(weights.size());
            r_sendbuff[owner].insert(r_sendbuff[owner].end(), normals.begin(),
                                     normals.end());
            r_sendbuff[owner].insert(r_sendbuff[owner].end(), weights.begin(),
                                     weights.end());
          }
          else
          {
            uint const glb_id = mesh.distdata().get_vertex_global(loc_id);
            dolfin_assert(glb_id < mesh.global_numVertices());
            VertexData * data = new VertexData();
            // Do not fill to avoid copy
            // getFacetData(weighting_, mesh, boundary, *bvertex,
            // data->facet_normals, data->facet_weights);
            vdata.insert(std::pair<uint, VertexData *>(glb_id, data));
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
    _set<uint> const& adjs = mesh.distdata().get_adj_ranks(0);
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      u_maxsendcount = std::max(u_maxsendcount, int(u_sendbuff[*it].size()));
      r_maxsendcount = std::max(r_maxsendcount, int(r_sendbuff[*it].size()));
    }
    dolfin_assert(u_maxsendcount <= u_size * mesh.topology().num_ghosts(0));
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
        dolfin_assert(glb_id < mesh.global_numVertices());
        uint const num_nc = u_recvbuff[iiu + 1];

        VertexDataMap::iterator it = vdata.find(glb_id);
        if (it == vdata.end())
        {
          error("Invalid adjacency data: unknown global vertex %d.", glb_id);
        }
        // Add corresponding facet normals and weights
        dolfin_assert(it->second != NULL);
        VertexData * vd = it->second;
        vd->facet_normals.insert(vd->facet_normals.end(), rptr,
                                 rptr + gdim * num_nc);
        rptr += gdim * num_nc;
        vd->facet_weights.insert(vd->facet_weights.end(), rptr, rptr + num_nc);
        rptr += num_nc;
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
  if (boundary.numCells() > 0)
  {
    boundary.init(boundary.topology().dim(), 0);
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
      uint const loc_id = boundary.vertex_index(*bvertex);
      if (mesh.distdata().is_ghost(loc_id, 0))
      {
        continue;
      }

      //--- Get facet normals and weights -----------------------------------
      Array<real> N;
      Array<real> W;
      getFacetData(type_, mesh, boundary, *bvertex, N, W);
      bool const is_shared = mesh.distdata().is_shared(loc_id, 0);
      uint const glb_id = mesh.distdata().get_vertex_global(loc_id);
      if (is_shared)
      {
        VertexDataMap::iterator it = vdata.find(glb_id);
        N.insert(N.end(), it->second->facet_normals.begin(),
                 it->second->facet_normals.end());
        W.insert(W.end(), it->second->facet_weights.begin(),
                 it->second->facet_weights.end());
      }
      dolfin_assert(N.size() == W.size() * gdim);

      //--- Compute basis ---------------------------------------------------
      uint num_surfaces = computeBasis(gdim, B, N, W, cosalpha, weighted);
      uint vertex_type = std::min(tdim, num_surfaces);

      //
      vertex_type_.set(loc_id, vertex_type);
      for (uint e = 0; e < gdim; ++e)
      {
        for (uint d = 0; d < gdim; ++d)
        {
          basis_[e][d].set(loc_id, B[e][d]);
        }
      }

      if (is_shared)
      {
        _set<uint> const& adjs = mesh.distdata().get_shared_adj(loc_id, 0);
        for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end();
             ++it)
        {
          u_sendbuff[*it].push_back(glb_id);
          u_sendbuff[*it].push_back(vertex_type);
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
    _set<uint> const& adjs = mesh.distdata().get_adj_ranks(0);
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
        dolfin_assert(mesh.distdata().has_global(u_recvbuff[iiu], 0));
        uint loc_id = mesh.distdata().get_vertex_local(u_recvbuff[iiu]);
        dolfin_assert(mesh.distdata().is_ghost(loc_id, 0));
        vertex_type_.set(loc_id, u_recvbuff[iiu + 1]);
        for (uint e = 0; e < gdim; ++e)
        {
          for (uint d = 0; d < gdim; ++d)
          {
            basis_[e][d].set(loc_id, r_recvbuff[iir]);
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
  for (VertexDataMap::iterator it = vdata.begin(); it != vdata.end(); ++it)
  {
    delete it->second;
  }
  vdata.clear();

}
//-----------------------------------------------------------------------------
uint VertexNormal::computeBasis(uint gdim, Point B[], Array<real> N,
                                Array<real> W, real cosalpha_max, bool weighted)
{
  if(N.size() != W.size() * gdim)
  {
    error("Mismatch between normals and weights: %d components, %d weights",
          N.size(), W.size());
  }

  //--- Determine vertex type by discriminating surfaces ----------------
  Array<Point> nS;
  Array<real> wS;
  std::set<uint> Rnormals;
  uint const num_facets = W.size();
  std::set<uint>::iterator it = Rnormals.begin();
  real wSa = 0.0;
  for (uint i = 0; i < num_facets; ++i)
  {
    Rnormals.insert(it, i);
    wSa += W[i];
    dolfin_assert(W[i] > 0.0);
  }

  while (Rnormals.size() > 0)
  {
    // Initialize new surface normal and weight with reference
    it = Rnormals.begin();
    uint const rfacet = (*it);
    Point nSx;
    real wSx = W[rfacet];
    for (uint d = 0; d < gdim; ++d)
    {
      nSx[d] = wSx * N[gdim * rfacet + d];
    }

    std::set<uint> Unormals;
    Unormals.insert(Unormals.begin(), rfacet);
    // Loop through remaining normals indexes
    for (++it; it != Rnormals.end(); ++it)
    {
      uint const cfacet = (*it);
      dolfin_assert(W[cfacet] > 0);

      // Compute the scalar product with the reference normal
      real cosalpha = 0.0;
#pragma _CRI novector
      for (uint d = 0; d < gdim; ++d)
      {
        cosalpha += N[gdim * rfacet + d] * N[gdim * cfacet + d];
      }
      if (cosalpha > cosalpha_max)
      {
        // Add contribution to surface normal
        for (uint d = 0; d < gdim; ++d)
        {
          nSx[d] += W[cfacet] * N[gdim * cfacet + d];
        }
        wSx += W[cfacet];
        // Take into account that the normal is used
        Unormals.insert(Unormals.begin(), cfacet);
      }
    }

    // Add unit surface normal to the list of surface normals
    //message("Surface normal %i", vertex_type);
    nSx /= nSx.norm();
    nS.push_back(nSx);
    if(weighted)
    {
      wSx /= wSa;
    }
    else
    {
      wSx = 1.0;
    }
    wS.push_back(wSx);

    // Next loop we add a new surface and discriminate again against
    // the remaining normals
    for (it = Unormals.begin(); it != Unormals.end(); ++it)
    {
      Rnormals.erase(*it);
    }
  }
  dolfin_assert(wS.size() > 0);

  //--- Compute vertex normal -------------------------------------------
  // The typical algorithm would be:
  // n_k    = sum_{i=1}^k nS_i
  // tau1_k = |n_{k-1}|^2 nS_k - (n_{k-1} . nS_k ) n_{k-1}
  // tau2_k = n_k ^ tau1_k
  // and such that in 2d tau2_k = ez = (0 , 0, 1)
  for (uint d = 0; d < gdim; ++d)
  {
    B[0][d] = 0.0;
    for (uint s = 0; s < nS.size(); ++s)
    {
      B[0][d] += wS[s] * nS[s][d];
    }
  }
  B[0] /= B[0].norm();

  switch (gdim)
    {
    case 1:
      break;
    case 2:
      B[1][0] = -B[0][1];
      B[1][1] = +B[0][0];
      break;
    case 3:
      if (nS.size() == 1)
      {
        if (std::fabs(B[0][0]) >= 0.5 || std::fabs(B[0][1]) >= 0.5)
        {
          // t1 = rotation in (x, y)
          B[1][0] = -B[0][1];
          B[1][1] = +B[0][0];
          B[1][2] = +0.0;
          B[1] /= B[1].norm();

          // t2 = n ^ t1
          B[2][0] = -B[0][2] * B[1][1];
          B[2][1] = +B[0][2] * B[1][0];
          B[2][2] = +B[0][0] * B[1][1] - B[0][1] * B[1][0];
          // B[2] is de facto normalized
        }
        else
        {
          // t1 = rotation in (y, z)
          B[1][0] = +0.0;
          B[1][1] = -B[0][2];
          B[1][2] = +B[0][1];
          B[1] /= B[1].norm();

          // t2 = n ^ t1
          B[2][0] = +B[0][1] * B[1][2] - B[0][2] * B[1][1];
          B[2][1] = -B[0][0] * B[1][2];
          B[2][2] = +B[0][0] * B[1][1];
          // B[2] is de facto normalized
        }
      }
      else
      {
        uint const k = nS.size() - 1;
        // t2 = n ^ nSk / || n ^ nSk ||
        B[2][0] = +B[0][1] * nS[k][2] - nS[k][1] * B[0][2];
        B[2][1] = +B[0][2] * nS[k][0] - nS[k][2] * B[0][0];
        B[2][2] = +B[0][0] * nS[k][1] - nS[k][0] * B[0][1];
        B[2] /= B[2].norm();

        // t1 = t2 ^ n
        B[1][0] = +B[2][1] * B[0][2] - B[2][2] * B[0][1];
        B[1][1] = +B[2][2] * B[0][0] - B[2][0] * B[0][2];
        B[1][2] = +B[2][0] * B[0][1] - B[2][1] * B[0][0];
        // B[1] is de facto normalized
      }
      break;
    default:
      error("Unsupported geometric dimension.");
      break;
    }

#if DEBUG
  for (uint i = 0; i < gdim; ++i)
  {
    for (uint d = 0; d < gdim; ++d)
    {
      if(B[i][d] != B[i][d])
      {
        error("Component %d of vector %d is Not-a-Number.", d, i);
      }
    }
    real en = B[i].norm();
    if(!abscmp(en, 1.0, gdim*DOLFIN_EPS))
    {
      error("Basis is not normal: ||e%u|| = %e", i, en);
    }
    uint j = (i + 1) % gdim;
    real sp = B[i].dot(B[j]);
    if(!abscmp(sp, 0.0, gdim*DOLFIN_EPS))
    {
      error("Basis is not orthogonal: e%u . e%u = %e", i, j, sp);
    }
  }
#endif

  return nS.size();
}

//-----------------------------------------------------------------------------

}

