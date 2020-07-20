// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/VertexNormal.h>

#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellIterator.h>
#include <dolfin/mesh/EuclideanBasis.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(VertexNormal& other) :
    mesh_(other.mesh_),
    gdim_(other.gdim_),
    subdomain_(nullptr),
    basis_(gdim_ * gdim_, MeshValues<real, Vertex>(mesh_)),
    vertex_type_(mesh_),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(none)
{
  vertex_type_ = other.vertex_type_;
  basis_ = other.basis_;
}

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(Mesh& mesh, Type weight) :
    mesh_(mesh),
    gdim_(mesh.geometry_dimension()),
    subdomain_(nullptr),
    basis_(gdim_ * gdim_, MeshValues<real, Vertex>(mesh_)),
    vertex_type_(mesh),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(weight)
{
  computeNormal(mesh);
}

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(Mesh& mesh, SubDomain const& subdomain, Type weight) :
    mesh_(mesh),
    gdim_(mesh.geometry_dimension()),
    subdomain_(&subdomain),
    basis_(gdim_ * gdim_, MeshValues<real, Vertex>(mesh_)),
    vertex_type_(mesh),
    alpha_max_(0.5 * DOLFIN_PI),
    type_(weight)
{
  computeNormal(mesh);
}

//-----------------------------------------------------------------------------
VertexNormal& VertexNormal::operator=(VertexNormal&)
{
  return *this;
}

//-----------------------------------------------------------------------------
void VertexNormal::getFacetData(VertexNormal::Type type, Mesh& mesh,
                                BoundaryMesh& boundary, Vertex& bvertex,
                                Array<real>& normals, Array<real>& weights)
{
  uint const tdim = mesh.topology_dimension();
  uint const gdim = mesh.geometry_dimension();
  for (CellIterator bcell(bvertex); !bcell.end(); ++bcell)
  {
    Facet facet(mesh, boundary.facet_index(*bcell));

    if(subdomain_ != nullptr && !subdomain_->inside(&bcell->midpoint()[0], true))
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

  uint const tdim = mesh.topology_dimension();
  uint const gdim = mesh.geometry_dimension();
  // Important: make sure facet to cell connectivities are initialized
  BoundaryMesh& boundary = mesh.exterior_boundary();

  VertexDataMap vdmap;
#ifdef HAVE_MPI
  int rank = dolfin::MPI::rank();
#endif
  int pe_size = dolfin::MPI::size();
  Array< Array<uint> > u_sendbuff( pe_size );
  Array< Array<real> > r_sendbuff( pe_size );

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
    int maxsendcount[2] = {0, 0};
    int maxrecvcount[2] = {0, 0};

    _set<uint> const & adjs = ddv.get_adj_ranks();
    for ( _set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it )
    {
      maxsendcount[0] = std::max(maxsendcount[0], int(u_sendbuff[*it].size()));
      maxsendcount[1] = std::max(maxsendcount[1], int(r_sendbuff[*it].size()));
    }
    dolfin_assert( ( uint ) maxsendcount[0] <= 2
                   and mesh.topology().num_ghost( 0 ) );

    MPI::all_reduce< MPI::max >( maxsendcount, maxrecvcount, 2 );

    dolfin_assert( maxrecvcount[0] > 0 );
    dolfin_assert( maxrecvcount[1] > 0 );

    // For each process
    Array< uint > u_recvbuff( maxrecvcount[0] );
    Array< real > r_recvbuff( maxrecvcount[1] );
    for ( int j = 1; j < pe_size; ++j )
    {
      uint src  = ( rank - j + pe_size ) % pe_size;
      uint dest = ( rank + j ) % pe_size;

      int u_recvcount = MPI::sendrecv( u_sendbuff[dest], dest,
                                       u_recvbuff, src, 1 );

      MPI::sendrecv( r_sendbuff[dest], dest, r_recvbuff, src, 1 );

      real * rptr = &r_recvbuff[0];
      for ( int iiu = 0; iiu < u_recvcount; iiu += 2 )
      {
        uint const glb_id = u_recvbuff[iiu];
        dolfin_assert( glb_id < mesh.global_size( 0 ) );
        uint const num_nc = u_recvbuff[iiu + 1];

        VertexDataMap::iterator it = vdmap.find( glb_id );
        if ( it != vdmap.end() )
        {
          // Add corresponding facet normals and weights
          dolfin_assert( it->second != nullptr );
          VertexData * vd = it->second;
          vd->facet_normals.insert( vd->facet_normals.end(), rptr,
                                    rptr + gdim * num_nc );
          rptr += gdim * num_nc;
          vd->facet_weights.insert( vd->facet_weights.end(), rptr,
                                    rptr + num_nc );
          rptr += num_nc;
        }
      }

      // Clear for reuse
      u_sendbuff[dest].clear();
      r_sendbuff[dest].clear();
    }
#endif
  }

  //--- Compute normals -------------------------------------------------------
  real const cosalpha = std::cos(alpha_max_);
  if (boundary.num_cells() > 0)
  {
    // Initialize cartesian basis
    Point B[Space::MAX_DIMENSION];
    for (uint d = 0; d < Space::MAX_DIMENSION; ++d)
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
      vertex_type_(local_index) = vtype;
      for (uint e = 0; e < gdim; ++e)
      {
        for (uint d = 0; d < gdim; ++d)
        {
          basis(e, d)(local_index) = B[e][d];
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

  if ( mesh.is_distributed() )
  {
#ifdef HAVE_MPI

    DistributedData const & ddv = mesh.distdata()[0];

    // Exchange data
    int maxsendcount[2] = {0, 0};
    int maxrecvcount[2] = {0, 0};

    _set<uint> const & adjs = ddv.get_adj_ranks();
    for ( _set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it )
    {
      maxsendcount[0] = std::max(maxsendcount[0], int(u_sendbuff[*it].size()));
      maxsendcount[1] = std::max(maxsendcount[1], int(r_sendbuff[*it].size()));
    }

    MPI::all_reduce< MPI::max >( maxsendcount, maxrecvcount, 2 );

    dolfin_assert( maxrecvcount[0] > 0 );
    dolfin_assert( maxrecvcount[1] > 0 );

    // For each process
    Array< uint > u_recvbuff( maxrecvcount[0] );
    Array< real > r_recvbuff( maxrecvcount[1] );

    for ( int j = 1; j < pe_size; ++j )
    {
      uint src  = ( rank - j + pe_size ) % pe_size;
      uint dest = ( rank + j ) % pe_size;

      int u_recvcount = MPI::sendrecv( u_sendbuff[dest], dest,
                                       u_recvbuff, src, 1 );

      MPI::sendrecv( r_sendbuff[dest], dest, r_recvbuff, src, 1 );

      //
      uint iir = 0;
      for ( int iiu = 0; iiu < u_recvcount; iiu += 2 )
      {
        dolfin_assert( ddv.has_global( u_recvbuff[iiu] ) );
        uint const local_index = ddv.get_local( u_recvbuff[iiu] );
        dolfin_assert( ddv.is_ghost( local_index ) );
        vertex_type_( local_index ) = u_recvbuff[iiu + 1];
        for ( uint e = 0; e < gdim; ++e )
        {
          for ( uint d = 0; d < gdim; ++d )
          {
            basis( e, d )( local_index ) = r_recvbuff[iir];
            ++iir;
          }
        }
      }
    }
#endif
  }

  // Cleanup
  for (VertexDataMap::iterator it = vdmap.begin(); it != vdmap.end(); ++it)
  {
    delete it->second;
  }
}

//-----------------------------------------------------------------------------

}
