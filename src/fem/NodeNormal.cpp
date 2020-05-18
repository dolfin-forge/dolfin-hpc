// Copyright (C) 2014-2016 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/NodeNormal.h>

#include <dolfin/common/timing.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/EuclideanBasis.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexNormal.h>

namespace dolfin
{

struct FacetData
{
  uint global_index;
  real weight;
  Point normal;
  _set<uint> nodes;

  void disp() const
  {
    section("FacetData");
    prm("global_index", global_index);
    prm("nodes"       , nodes.size());
    prm("weight"      , weight);
    prm("normal"      , normal);
    end();
  }
};

struct NodeData
{
  uint node_type;
  Array<uint> dofs;
  Array<uint> adjs;
  Array<FacetData *> facets;

  void disp() const
  {
    section("NodeData");
    prm("node_type" , node_type);
    prm("dofs"      , dofs.size());
    prm("adjs"      , adjs.size());
    prm("facets"    , facets.size());
    end();
  }
};

typedef _map<uint, FacetData *>::iterator FacetMapIterator;
typedef _map<uint, NodeData *>::iterator  NodeMapIterator;

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal( Mesh & mesh, Type w, real alpha )
  : BoundaryNormal( mesh )
  , mesh_( mesh )
  , subdomain_( NULL )
  , alpha_max_( alpha )
  , type_( w )
{
}

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal( Mesh & mesh,
                        SubDomain const & subdomain,
                        Type w,
                        real alpha )
  : BoundaryNormal( mesh )
  , mesh_( mesh )
  , subdomain_( &subdomain )
  , alpha_max_( alpha )
  , type_( w )
{
}

//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
}

//-----------------------------------------------------------------------------
void NodeNormal::clear()
{
  node_type_.clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::compute()
{
  uint const gdim = mesh_.geometry_dimension();
  if (basis().size() < gdim)
  {
    error("Invalid size of storage vector for basis functions in NodeNormal");
  }
  compute(mesh_, basis());
}

//-----------------------------------------------------------------------------
uint NodeNormal::node_type(uint node_id) const
{
  dolfin_assert(node_type_.size() > 0);
  _map<uint, uint>::const_iterator it = node_type_.find(node_id);
  if (it == node_type_.end())
  {
    error("Invalid node id requested for node type");
  }
  return it->second;
}

//-----------------------------------------------------------------------------
void NodeNormal::compute(Mesh& mesh, Array<Function>& basis)
{
  message(1, "NodeNormal : compute P%u node normal", basis[0].space().degree());
  clear();
  BoundaryMesh& boundary = mesh.exterior_boundary();

  tic();
  //---------------------------------------------------------------------------
  uint const tdim = mesh.topology_dimension();
  uint const fdim = mesh.type().facet_dim();
  uint const gdim = mesh.geometry_dimension();
#ifdef HAVE_MPI
  uint const pe_size = MPI::size();
  uint const rank    = MPI::rank();
#endif

  // Maps facet global index to (weight, normal)
  _map<uint, FacetData *> facets_data;
  facets_data.reserve( boundary.num_vertices() );
  // Maps dofs to facet global indices
  _map<uint, NodeData *> nodes_data;
  nodes_data.reserve( boundary.num_cells() );

  //[facet, nb_nodes, [node indices]]
  Array< Array<uint> > u_sendbuf( pe_size );
  //[weight, normal]
  Array< Array<real> > r_sendbuf( pe_size );

  FiniteElementSpace const& spaceN = basis[0].space();
  DofMap const& dofmapN = spaceN.dofmap();
  ScratchSpace scratchN(spaceN);
  uint const value_size = spaceN.element().value_size();
  dolfin_assert(value_size == gdim);
  uint const num_facet_dofs = dofmapN.num_facet_dofs();
  uint const num_facet_nodes = num_facet_dofs / value_size;
  uint const num_restricted_facet_dofs = dofmapN.num_entity_dofs(fdim);

  // The implementation works only for dofs located on the exterior boundary
  bool const on_boundary = true;

  //--- Create facet data
  // Mark facets in the subdomain based on dofs, naive implementation
  node_type_.clear();
  _set<uint> used_ghost_nodes;
  for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
  {
    Facet facet(mesh, boundary.facet_index(*bcell));
    Cell cell(mesh, facet.entities(tdim)[0]);
    uint local_facet = cell.index(facet);

    // An exterior facet should be included in the subdomain if at least one
    // of the dofs on the facet restriction (if any) or if the facet midpoint
    // is in the subdomain.
    // Skip the facet if it does not satisfy one of these conditions.
    if ((subdomain_ == NULL)
        || subdomain_->inside(&(bcell->midpoint())[0], on_boundary))
    {
      scratchN.cell.update(cell);
      dofmapN.tabulate_coordinates(scratchN.coordinates, scratchN.cell);
    }
    else
    {
      bool invalid = true;
      if (num_restricted_facet_dofs > 0)
      {
        // Tabulate dofs on facet restriction
        scratchN.cell.update(cell);
        dofmapN.tabulate_coordinates(scratchN.coordinates, scratchN.cell);
        dofmapN.tabulate_entity_dofs(scratchN.facet_dofs, fdim, local_facet);
        for (uint i = 0; i < num_restricted_facet_dofs; ++i)
        {
          uint loc_dof = scratchN.facet_dofs[i];
          if (subdomain_->inside(scratchN.coordinates[loc_dof], on_boundary))
          {
            invalid = false;
            break;
          }
        }
      }
      if (invalid)
      {
        // Skip facet
        continue;
      }
    }

    // Add facet to the list with facet weight and normal
    FacetData * data = new FacetData();
    data->global_index = facet.global_index();
    switch (type_)
      {
      case NodeNormal::facet:  // facet area
        data->weight = bcell->volume();
        break;
      case NodeNormal::none:  // no weight
      case NodeNormal::unit:  // unit per facet
      default:
        data->weight = 1.0;
        break;
      }
    data->normal = cell.normal(local_facet);
    dolfin_assert( abscmp( data->normal.norm(), 1.0, 1e-13 ) );

    // Set valid dofs within the current facet
    _ordered_set<uint> ghost_nodes;
    dofmapN.tabulate_dofs(scratchN.dofs, scratchN.cell, cell);
    dofmapN.tabulate_facet_dofs(scratchN.facet_dofs, local_facet);
    for (uint f_n = 0; f_n < num_facet_nodes; ++f_n)
    {
      uint dof0 = scratchN.facet_dofs[f_n];
      if ((subdomain_ == NULL)
          || subdomain_->inside(scratchN.coordinates[dof0], on_boundary))
      {
        // Take global dof index of the first component as node id
        uint const node_id = scratchN.dofs[dof0];
        data->nodes.insert(node_id);
        node_type_[node_id] = 0;
        if (!dofmapN.is_ghost(node_id))
        {
          // Trigger creation of owned node data if the entry does not exist
          NodeMapIterator it = nodes_data.find(node_id);
          if (it == nodes_data.end())
          {
            NodeData * n_data = new NodeData();
            for (uint d = 0; d < value_size; ++d)
            {
              uint local_dof = scratchN.facet_dofs[d * num_facet_nodes + f_n];
              n_data->dofs.push_back(scratchN.dofs[local_dof]);
            }
            n_data->facets.push_back(data);
            nodes_data[node_id] = n_data;
          }
          else
          {
            it->second->facets.push_back(data);
          }
        }
        // Add to the set of ghosted dofs for sending data
        else if(used_ghost_nodes.count(node_id) == 0)
        {
          ghost_nodes.insert(node_id);
        }
      }
    }
    dolfin_assert(!data->nodes.empty());
    dolfin_assert(data->nodes.size() == num_facet_nodes);
    facets_data.insert(std::pair<uint, FacetData *>(data->global_index, data));

    // Collect to send ghosted data to the adjacent ranks
    /// @todo Dirty hack, send to all adjacents, not optimal but good enough
    if (!ghost_nodes.empty())
    {
      _set<uint> adjs;
      adjs.reserve( pe_size );
      for (VertexIterator v(facet); !v.end(); ++v)
      {
        if(v->is_shared())
        {
          _set<uint> const& a = mesh.distdata()[0].get_shared_adj(v->index());
          adjs.insert(a.begin(), a.end());
        }
      }
      for(_set<uint>::const_iterator ai = adjs.begin(); ai != adjs.end(); ++ai)
      {
        // dofs
        dolfin_assert(ghost_nodes.size() <= num_facet_nodes);
        u_sendbuf[*ai].push_back(ghost_nodes.size());
        u_sendbuf[*ai].append(ghost_nodes.begin(), ghost_nodes.end());
        // global index
        u_sendbuf[*ai].push_back(data->global_index);
        // weight
        r_sendbuf[*ai].push_back(data->weight);
        // normal
        for (uint d = 0; d < gdim; ++d)
        {
          dolfin_assert(std::fabs(data->normal[d]) < 1.0 + DOLFIN_EPS);
          r_sendbuf[*ai].push_back(data->normal[d]);
        }
      }
      used_ghost_nodes.insert(ghost_nodes.begin(), ghost_nodes.end());
    }
  }

  //--- Exchange data for exterior facets with shared entities
  if (mesh.is_distributed())
  {
#ifdef HAVE_MPI
    // Since an entity is shared if all it lower dimensional entities are shared
    // we can loop over shared vertices and stack facets. If non-matching facet
    // are send they will be eventually discarded. This does not hold if the
    // subdomain has a hole in the interior of the facet.

    // Exchange values
    int maxsendcount[2] = { 0, 0 };
    int maxrecvcount[2] = { 0, 0 };

    _set<uint> const& adjs = mesh.distdata()[0].get_adj_ranks();
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      maxsendcount[0] = std::max( maxsendcount[0], int(u_sendbuf[*it].size()) );
      maxsendcount[1] = std::max( maxsendcount[1], int(r_sendbuf[*it].size()) );
    }

    MPI::all_reduce<MPI::max>( maxsendcount, maxrecvcount, 2 );

    dolfin_assert( maxrecvcount[0] >= maxsendcount[0] );
    dolfin_assert( maxrecvcount[1] >= maxsendcount[1] );

    Array< uint > u_recvbuf( maxrecvcount[0] );
    Array< real > r_recvbuf( maxrecvcount[1] );

    for (uint j = 1; j < pe_size; ++j)
    {
      int src  = (rank - j + pe_size) % pe_size;
      int dest = (rank + j) % pe_size;

      int n_u_recv = MPI::sendrecv( u_sendbuf[dest], dest, u_recvbuf, src, 1);
      MPI::sendrecv( r_sendbuf[dest], dest, r_recvbuf, src, 1 );

      uint iir = 0;
      for ( int iiu = 0; iiu < n_u_recv; )
      {
        Array<NodeMapIterator> valid;
        uint nb_nodes = u_recvbuf[iiu++];
        dolfin_assert(nb_nodes <= num_facet_nodes);
        // Check if one node is owned
        for (uint n = 0; n < nb_nodes; ++n)
        {
          NodeMapIterator it = nodes_data.find(u_recvbuf[iiu++]);
          if (it != nodes_data.end())
          {
            // Process owns the node
            valid.push_back(it);
            dolfin_assert(!dofmapN.is_ghost(it->first));
          }
        }
        uint global_index = u_recvbuf[iiu++];

        if(!valid.empty())
        {
          dolfin_assert(mesh.distdata()[0].get_adj_ranks().count(src) > 0);
          // Add facet to adjacent process
          FacetData * data = new FacetData();
          data->global_index = global_index;
          data->weight = r_recvbuf[iir];
          std::copy(&r_recvbuf[iir+1], &r_recvbuf[iir+1] + gdim, &data->normal[0]);

          // Add facet to the nodes' list of adjacent facets
          for(uint n = 0; n < valid.size(); ++n)
          {
            data->nodes.insert(valid[n]->first);
            valid[n]->second->facets.push_back(data);
            valid[n]->second->adjs.push_back(src);
          }

          // Add facet
          facets_data.insert(std::pair<uint, FacetData *>(data->global_index, data));
        }
        iir += gdim + 1;
        dolfin_assert( iiu <= n_u_recv );
      }

      // Clear for reuse
      u_sendbuf[dest].clear();
      r_sendbuf[dest].clear();
    }
#endif // HAVE_MPI
  }

  //--- Determine node type from facet normals and compute surface normals
  real const cosalpha = std::cos(alpha_max_);
  bool const weighted = (type_ != NodeNormal::none);
  uint const num_boundary_dofs = gdim * nodes_data.size();

  Array< uint > dofs( num_boundary_dofs );  // dof indices
  uint node_dofs = 0;

  Array< real > block( gdim * num_boundary_dofs );  // ( n, tau_1, tau_2 )
  uint offset = 0;

  // Initialize cartesian basis
  Point B[Space::MAX_DIMENSION];

  Function nn;
  Array< real > nnblock;

  if ( dolfin_get< bool >( "io_print_node_normal_types" ) == true )
  {
    nn = basis[0];
    nn = 100;
    nnblock.resize( num_boundary_dofs );
  }

  for (uint d = 0; d < Space::MAX_DIMENSION; ++d)
  {
    B[d][d] = 1.0;
  }
  // Contains only owned nodes
  for ( NodeMapIterator it = nodes_data.begin();
        it != nodes_data.end();
        ++it, node_dofs += gdim, offset+=gdim )
  {
    uint const node_id = it->first;
    NodeData * n_data = it->second;
    dolfin_assert(!dofmapN.is_ghost(node_id));

    // Collect facets data
    Array<real> N;
    Array<real> W;
    for(uint i = 0; i < n_data->facets.size(); ++i)
    {
      FacetData * f_data = n_data->facets[i];
      // Ghosted facets data only contains ghosted nodes
      dolfin_assert(f_data->nodes.size() <= num_facet_nodes);
      N.insert(N.end(), &f_data->normal[0], &f_data->normal[0] + gdim);
      W.push_back(f_data->weight);
    }

    // Compute basis and node type
    uint numS = EuclideanBasis::compute(gdim, B, N, W, cosalpha, weighted);
    uint node_type = std::min(tdim, numS);
    if ( dolfin_get<bool>( "NodeNormal restricted" ) == true )
      node_type = ( node_type == 2 ) ? 3 : node_type;

    // Set node type
    node_type_[node_id] = node_type;
    n_data->node_type = node_type;// Useless if node data is not reused
    dolfin_assert(n_data->node_type > 0);
    // If the node is shared, prepare for synchronization with adjacents
    for(Array<uint>::const_iterator ait = n_data->adjs.begin();
    ait != n_data->adjs.end(); ++ait)
    {
      u_sendbuf[*ait].push_back(node_id);
      u_sendbuf[*ait].push_back(node_type);
    }

    // Copy dof indices to array for vector block set.
    std::copy(n_data->dofs.begin(),n_data->dofs.end(), dofs.data() + node_dofs );

    if ( dolfin_get< bool >( "io_print_node_normal_types" ) == true )
    {
      nnblock[offset+0] = node_type;
      nnblock[offset+1] = node_id;
      nnblock[offset+2] = 1.0 * dofmapN.is_shared(node_id);
    }

    // Copy data to block array
    for(uint d = 0; d< gdim; ++d)
    {
      std::copy(&B[d][0], &B[d][0] + gdim, block.data() + offset + d*num_boundary_dofs);
    }
  }

  if ( dolfin_get< bool >( "io_print_node_normal_types" ) == true )
  {
    GenericVector& v = nn.vector();
    v.set( nnblock.data(), num_boundary_dofs, dofs.data() );
    v.apply();
    File( "node_normal_types.pvd" ) << nn;
  }

  // Set vectors, values are synchronized
  for (uint d = 0; d < gdim; ++d)
  {
    GenericVector& v = basis[d].vector();
    v.set( block.data() + d * num_boundary_dofs, num_boundary_dofs, dofs.data() );
    v.apply();
  }

  // Cleanup facets and nodes data
  for (FacetMapIterator it = facets_data.begin();
  it != facets_data.end(); ++it)
  {
    delete it->second;
  }
  facets_data.clear();
  for (NodeMapIterator it = nodes_data.begin();
  it != nodes_data.end(); ++it)
  {
    delete it->second;
  }
  nodes_data.clear();

  // Synchronize node type
  if (mesh.is_distributed())
  {
#if HAVE_MPI
    uint const u_size = 2;
    int u_maxsendcount = 0;
    int u_maxrecvcount = 0;
    _set<uint> const& adjs = mesh.distdata()[0].get_adj_ranks();
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      u_maxsendcount = std::max(u_maxsendcount, int(u_sendbuf[*it].size()));
    }
    MPI::all_reduce<MPI::max>( u_maxsendcount, u_maxrecvcount );
    dolfin_assert(u_maxrecvcount > 0);

    // For each process
    Array< uint > u_recvbuf( u_maxrecvcount );
    for (uint j = 1; j < pe_size; ++j)
    {
      int src = (rank - j + pe_size) % pe_size;
      int dst = (rank + j) % pe_size;

      int n_u_recv = MPI::sendrecv( u_sendbuf[dst], dst, u_recvbuf, src, 1 );

      // Add node type to the map
      for (int iiu = 0; iiu < n_u_recv; iiu += u_size)
      {
        uint const id = u_recvbuf[iiu];
        uint const nt = u_recvbuf[iiu + 1];

        dolfin_assert(node_type_[id] == 0); // Should be sent once
        dolfin_assert(dofmapN.is_ghost(id));
        dolfin_assert(nt > 0);

        node_type_[id] = nt;
      }
    }
#endif
  }

  tocd(1);

  for (uint e = 0; e < gdim; ++e)
  {
    basis[e].sync();
  }
}

//-----------------------------------------------------------------------------

}
