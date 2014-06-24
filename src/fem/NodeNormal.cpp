// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2009.
// Modified by Aurélien Larcher, 2012-14. (rewrite, extension to any element)
//
// First added:  2007-05-01
// Last changed: 2014-05-22

#include <dolfin/fem/NodeNormal.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>

#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Mesh& mesh, Type w, real alpha) :
    BoundaryNormal(mesh),
    mesh_(mesh),
    subdomain_(NULL),
    no_subdomain_(true),
    alpha_max_(alpha),
    weighting_(w)
{
}

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Mesh& mesh, SubDomain const& subdomain, Type w,
                       real alpha) :
    BoundaryNormal(mesh),
    mesh_(mesh),
    subdomain_(&subdomain),
    no_subdomain_(false),
    alpha_max_(alpha),
    weighting_(w)
{
}

//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
  Clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::Clear()
{
  // Cleanup facet and node data
  for (_map<uint, FacetData *>::iterator it = facets_.begin();
  it != facets_.end(); ++it)
  {
    delete it->second;
  }
  facets_.clear();
  for (_map<uint, NodeData *>::iterator it = nodes_.begin();
      it != nodes_.end(); ++it)
  {
    delete it->second;
  }
  nodes_.clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::compute()
{
  Compute(mesh_, basis());
}

//-----------------------------------------------------------------------------
uint NodeNormal::node_type(uint node_id) const
{
  dolfin_assert(nodes_.size() > 0);
  _map<uint, NodeData *>::const_iterator it = nodes_.find(node_id);
  dolfin_assert(it != nodes_.end());
  return it->second->node_type;
}

//-----------------------------------------------------------------------------
void NodeNormal::Compute(Mesh& mesh, Array<Function>& functions)
{
  Clear();

  BoundaryMesh& boundary = mesh.exterior_boundary();
  if (!boundary.numCells() > 0)
  {
    return;
  }
  MeshFunction<uint> * cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint> * vertex_map = boundary.data().meshFunction("vertex map");

  //---------------------------------------------------------------------------
  MeshDistributedData & distdata = mesh.distdata();
  uint const tdim = mesh.topology().dim();
  uint const facet_dim = tdim - 1;
  uint const gdim = mesh.geometry().dim();

  //
  if (functions.size() < gdim)
  {
    error("Invalid size of storage vector for basis functions in NodeNormal");
  }
  for (uint d = 0; d < gdim; ++d)
  {
    if (functions[d].type() != Function::discrete)
    {
      error("All basis functions in NodeNormal should be discrete");
    }
  }
  FiniteElementSpace const& space = functions[0].space();
  DofMap const& dofmap = space.dofmap();
  ScratchSpace scratch(space);
  uint const value_dim = space.element().value_dimension(0);
  dolfin_assert(value_dim == gdim);
  uint const num_facet_dofs = dofmap.num_facet_dofs();
  uint const num_facet_nodes = num_facet_dofs / value_dim;
  uint const num_restricted_facet_dofs = dofmap.num_entity_dofs(facet_dim);

  // The implementation works only for dofs located on the exterior boundary
  bool const on_boundary = true;

  //--- Create facet data
  // Mark facets in the subdomain based on dofs, naive implementation
  for (CellIterator b_cell(boundary); !b_cell.end(); ++b_cell)
  {
    Facet facet(mesh, cell_map->get(*b_cell));
    Cell cell(mesh, facet.entities(tdim)[0]);
    uint local_facet = cell.index(facet);

    // An exterior facet should be included in the subdomain if at least one
    // of the dofs on the facet restriction (if any) or if the facet midpoint
    // is in the subdomain.
    // Skip the facet if it does not satifies one of these conditions.
    if (no_subdomain_
        || subdomain_->inside(&(b_cell->midpoint())[0], on_boundary))
    {
      // Update cell data and tabulate the coordinates
      scratch.cell.update(cell, distdata);
      dofmap.tabulate_coordinates(scratch.coordinates, scratch.cell);
    }
    else
    {
      bool invalid = true;
      if (num_restricted_facet_dofs > 0)
      {
        // Tabulate dofs on facet restriction
        scratch.cell.update(cell, distdata);
        dofmap.tabulate_coordinates(scratch.coordinates, scratch.cell);
        dofmap.tabulate_entity_dofs(scratch.facet_dofs, facet_dim, local_facet);
        for (uint i = 0; i < num_restricted_facet_dofs; ++i)
        {
          uint loc_dof = scratch.facet_dofs[i];
          if (subdomain_->inside(scratch.coordinates[loc_dof], on_boundary))
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
    data->global_index = distdata.get_facet_global(facet.index());
    switch (weighting_)
    {
      case NodeNormal::none: // unit
        data->weight = 1.0;
        break;
      case NodeNormal::facet: // facet area
        data->weight = b_cell->volume();
        break;
      case NodeNormal::cell: // adjacent cell volume
        data->weight = cell.volume();
        break;
      default:
        break;
    }
    data->normal = cell.normal(local_facet);

    // Set valid dofs within the current facet
    dofmap.tabulate_dofs(scratch.dofs, scratch.cell, cell.index());
    dofmap.tabulate_facet_dofs(scratch.facet_dofs, local_facet);
    for (uint f_n = 0; f_n < num_facet_nodes; ++f_n)
    {
      uint dof0 = scratch.facet_dofs[f_n];
      if (no_subdomain_
          || subdomain_->inside(scratch.coordinates[dof0], on_boundary))
      {
        // Take global dof index of the first component as node id
        uint node_id = scratch.dofs[dof0];
        data->nodes.insert(node_id);

        // Trigger creation of node data if the entry does not exist
        _map<uint, NodeData *>::iterator it = nodes_.find(node_id);
        if (it == nodes_.end())
        {
          NodeData * n_data = new NodeData();
          for (uint d = 0; d < value_dim; ++d)
          {
            uint local_dof = scratch.facet_dofs[d * num_facet_nodes + f_n];
            n_data->dofs.push_back(scratch.dofs[local_dof]);
            dolfin_assert(no_subdomain_
                || subdomain_->inside(scratch.coordinates[local_dof],
                    on_boundary));
          }
          n_data->facets.push_back(data);
          nodes_[node_id] = n_data;
        }
        else
        {
          it->second->facets.push_back(data);
        }
      }
    } dolfin_assert(!data->nodes.empty());
    facets_.insert(std::pair<uint, FacetData *>(data->global_index, data));
  }

  //--- Exchange data for exterior facets with shared entities
  if (mesh.is_distributed())
  {
    // Since an entity is shared is shared iff all it lower dimensional entities
    // are shared we can loop over shared vertices and stack facets.
    // If non-matching facet are send they will be eventually discarded.
    // This does not hold if the subdomain has a hole in the interior of the
    // facet.
    uint rank = dolfin::MPI::processNumber();
    uint pe_size = dolfin::MPI::numProcesses();
    Array<uint> u_sendbuf;    //[facet, nb_nodes, [node indices]]
    Array<real> r_sendbuf;    //[weight, normal]
    uint const r_packet_size = 1 + gdim;
    _set<uint> used_adj_facets;
    for (VertexIterator boundary_vertex(boundary); !boundary_vertex.end();
        ++boundary_vertex)
    {
      uint vertex_idx = vertex_map->get(*boundary_vertex);
      if (!mesh.distdata().is_shared(vertex_idx, 0))
      {
        continue;
      }

      Vertex v(mesh, vertex_idx);
      uint const num_adj_facets = v.numEntities(facet_dim);
      dolfin_assert(num_adj_facets > 0);
      uint * adj_facets_idx = v.entities(facet_dim);
      for (uint f = 0; f < num_adj_facets; ++f)
      {
        uint const f_local = adj_facets_idx[f];

        // Avoid sending the same facet twice
        if (used_adj_facets.count(f_local) > 0)
        {
          continue;
        }
        used_adj_facets.insert(f_local);

        // Pack data to send buffer
        uint const f_global = distdata.get_facet_global(f_local);
        _map<uint, FacetData *>::iterator it = facets_.find(f_global);
        if (it != facets_.end())
        {
          // Append data as the facet shares a dof
          FacetData * data = it->second;
          // global index
          u_sendbuf.push_back(data->global_index);
          // dofs
          dolfin_assert(data->nodes.size() > 0);
          u_sendbuf.push_back(data->nodes.size());
          for (_set<uint>::const_iterator d_it = data->nodes.begin();
          d_it != data->nodes.end(); ++d_it)
          {
            u_sendbuf.push_back(*d_it);
          }
          // weight
          r_sendbuf.push_back(data->weight);
          // normal
          for (uint d = 0; d < gdim; ++d)
          {
            r_sendbuf.push_back(data->normal[d]);
          }
        }
      }
    }

    // Exchange values
    MPI_Status status;
    uint src;
    uint dest;
    MPI_Barrier(dolfin::MPI::DOLFIN_COMM);

    int  u_sendcount = u_sendbuf.size();
    int u_maxrecvcount = 0;
    int u_recvcount = 0;
    MPI_Allreduce(&u_sendcount, &u_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    uint * u_recvbuf = new uint[u_maxrecvcount];

    int r_sendcount = r_sendbuf.size();
    int r_maxrecvcount = 0;
    int r_recvcount = 0;
    MPI_Allreduce(&r_sendcount, &r_maxrecvcount, 1, MPI_INT, MPI_MAX,
                  dolfin::MPI::DOLFIN_COMM);
    real * r_recvbuf = new real[r_maxrecvcount];

    for (int proc = 1; proc < pe_size; ++proc)
    {
      src = (rank - proc + pe_size) % pe_size;
      dest = (rank + proc) % pe_size;

      MPI_Sendrecv(&u_sendbuf[0], u_sendcount, MPI_UNSIGNED, src, 1, u_recvbuf,
                   u_maxrecvcount, MPI_UNSIGNED, dest, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &u_recvcount);
      MPI_Sendrecv(&r_sendbuf[0], r_sendcount, MPI_DOUBLE, src, 1, r_recvbuf,
                   r_maxrecvcount, MPI_DOUBLE, dest, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &r_recvcount);

      uint i_r = 0;
      for (uint i_u = 0; i_u < u_recvcount;)
      {
        // Add adjacent process facet data
        FacetData * data = new FacetData();
        data->global_index = u_recvbuf[i_u++];
        uint nb_nodes = u_recvbuf[i_u++];
        dolfin_assert(nb_nodes > 0);dolfin_assert(nb_nodes <= num_facet_nodes);
        data->nodes.insert(&u_recvbuf[i_u], &u_recvbuf[i_u + nb_nodes]);
        i_u += nb_nodes;
        data->weight = r_recvbuf[i_r++];
        std::copy(&r_recvbuf[i_r], &r_recvbuf[i_r + gdim], &data->normal[0]);
        i_r += gdim;
        facets_.insert(std::pair<uint, FacetData *>(data->global_index, data));

        // Add facet to the nodes' list of adjacent facets
        for (_set<uint>::const_iterator d_it = data->nodes.begin();
        d_it != data->nodes.end(); ++d_it)
        {
          // If the node exists on the partition then add the facet
          _map<uint, NodeData *>::iterator it = nodes_.find(*d_it);
          if(it != nodes_.end())
          {
            it->second->facets.push_back(data);
          }
        }
      }
    }
    delete[] r_recvbuf;
    delete[] u_recvbuf;
  }

  //--- Determine node type from facet normals and compute surface normals
  real cosalpha_max = std::cos(alpha_max_);
  real cosalpha = 0.0;
  uint const num_boundary_dofs = gdim * nodes_.size();
  uint * dofs = new uint[num_boundary_dofs]; // dof indices
  uint * node_dofs = &dofs[0];
  real * block = new real[gdim * num_boundary_dofs]; // ( n, tau_1, tau_2 )
  real * offset = &block[0];
  Point basis[3]; //TODO: Remove Point
  for (_map<uint, NodeData *>::iterator it = nodes_.begin();
  it != nodes_.end(); ++it, node_dofs+=gdim, offset+=gdim)
  {
    NodeData * n_data = it->second;
    // Copy dof indices to array for vector block set.
    std::copy(n_data->dofs.begin(),n_data->dofs.end(),node_dofs);

    //
    Array<std::pair<real, Point> > surfaces;
    Array<FacetData *>const & n_facets = n_data->facets;
    uint const num_facets = n_facets.size();

    uint num_remaining_facets = num_facets;
    FacetData ** remaining_facets = new FacetData*[num_remaining_facets];
    std::copy(n_facets.begin(), n_facets.end(), remaining_facets);
    while(num_remaining_facets > 0)
    {
      // Set reference normal as the last of the remaining facet normal
      // and initialize new surface.
      Point& reference_normal = remaining_facets[0]->normal;
      real s_weight = remaining_facets[0]->weight;
      Point s_normal(reference_normal);
      uint num_eliminated_facets = 1;

      //
      uint entry_to_update = 0;
      for(uint f = 1; f < num_remaining_facets; ++f)
      {
        FacetData * curr_facet = remaining_facets[f];
        cosalpha = curr_facet->normal.dot(reference_normal);
        if (cosalpha > cosalpha_max)
        {
          // Add facet to current surface
          s_weight += curr_facet->weight;
          s_normal += curr_facet->weight * curr_facet->normal;
          ++num_eliminated_facets;
        }
        else
        {
          // Update current cursor with the facet for next loop
          remaining_facets[entry_to_update] = curr_facet;
          ++entry_to_update;
        }
      }
      s_weight /= num_eliminated_facets;
      s_normal /= s_normal.norm();
      surfaces.push_back(std::pair<real, Point>(s_weight, s_normal));
      num_remaining_facets -= num_eliminated_facets;

    }
    // Update node type with the number of discriminated surfaces
    n_data->node_type = surfaces.size();
    dolfin_assert(n_data->node_type > 0);

    //--- Compute node normals for piecewise linear boundary
    basis[0] = 0;
    for(Array<std::pair<real, Point> >::const_iterator s_it = surfaces.begin();
    s_it != surfaces.end(); ++s_it)
    {
      basis[0] += s_it->first * s_it->second;
    }
    basis[0] /= basis[0].norm();

    // Compute tangential vectors
    switch(gdim)
    {
      case 2:
      ComputeTangents2D(basis);
      break;
      case 3:
      if(n_data->node_type == 1)
      {
        ComputeTangents3DSurface(basis);
      }
      else
      {
        ComputeTangents3D(basis, surfaces.rbegin()->second);
      }
      break;
      default:
      break;
    }

    // Copy data to block array, FIXME: avoid copy by removing Points
    for(uint d = 0; d< gdim; ++d)
    {
      dolfin_assert(std::fabs(basis[d].norm() - 1.0) < DOLFIN_EPS);
      dolfin_assert(std::fabs(basis[d].dot(basis[(d+1)%(gdim)])) < DOLFIN_EPS);

      std::copy(&basis[d][0], &basis[d][0]+gdim, offset + d*num_boundary_dofs);
    }
  }

  for (uint d = 0; d < gdim; ++d)
  {
    GenericVector& v = functions[d].vector();
    v.set(block + d * num_boundary_dofs, num_boundary_dofs, dofs);
    v.apply();
  }
  delete[] dofs;
  delete[] block;
}

//-----------------------------------------------------------------------------
void NodeNormal::ComputeTangents2D(Point (&basis)[3])
{
  basis[1][0] = -basis[0][1];
  basis[1][1] = +basis[0][0];
}

//-----------------------------------------------------------------------------
void NodeNormal::ComputeTangents3DSurface(Point (&basis)[3])
{
  real norm_inv = 0.0;
  if (std::fabs(basis[0][0]) >= 0.5 || std::fabs(basis[0][1]) >= 0.5)
  {
    norm_inv = 1.
        / std::sqrt(basis[0][0] * basis[0][0] + basis[0][1] * basis[0][1]);
    // t11 = n2/n
    basis[1][0] = basis[0][1] * norm_inv;
    // t12 = -n1/n
    basis[1][1] = -basis[0][0] * norm_inv;
    // t13 = 0
    basis[1][2] = 0.0;
    // t21 = -t12*n3
    basis[2][0] = -basis[1][1] * basis[0][2];
    // t22 = t11*n3
    basis[2][1] = basis[1][0] * basis[0][2];
    // t23 = t12*n1 - t11*n2
    basis[2][2] = basis[1][1] * basis[0][0] - basis[1][0] * basis[0][1];
  }
  else
  {
    norm_inv = 1.
        / std::sqrt(basis[0][1] * basis[0][1] + basis[0][2] * basis[0][2]);
    // t11 = 0
    basis[1][0] = 0.0;
    // t12 = -n3/n
    basis[1][1] = -basis[0][2] * norm_inv;
    // t13 = n2/n
    basis[1][2] = basis[0][1] * norm_inv;
    // t21 = t13*n2 - t12*n3
    basis[2][0] = basis[1][2] * basis[0][1] - basis[1][1] * basis[0][2];
    // t22 = -t13*n1
    basis[2][1] = -basis[1][2] * basis[0][0];
    // t23 = t12*n1
    basis[2][2] = basis[1][1] * basis[0][0];
  }
}

//-----------------------------------------------------------------------------
void NodeNormal::ComputeTangents3D(Point (&basis)[3], Point& surface)
{
  basis[2] = basis[0].cross(surface);
  basis[2] /= basis[2].norm();
  basis[1] = basis[2].cross(basis[0]);
}
//-----------------------------------------------------------------------------

}

