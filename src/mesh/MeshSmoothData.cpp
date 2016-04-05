// Copyright (C) 2011 Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2011.
//
// First added:  2011-06-30
// Last changed: 2011-06-30

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/MeshSmoothData.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshSmoothData::MeshSmoothData(Mesh& mesh) :
    mesh_(mesh),
    boundary_(NULL)
{
}

//-----------------------------------------------------------------------------
MeshSmoothData::~MeshSmoothData()
{
}

//-----------------------------------------------------------------------------
void MeshSmoothData::prepare_mesh()
{
  if (boundary_ != NULL)
  {
    error("In MeshSmoothData, calling prepare_mesh() twice.");
  }

  //mapping for different boundaries

  //global boundary
  BoundaryMesh& boundary_global = mesh_.exterior_boundary();
  on_boundary_global_.init(mesh_, 0);
  for (VertexIterator vertex(boundary_global); !vertex.end(); ++vertex)
  {
    on_boundary_global_.set(boundary_global.vertex_index(*vertex), true);
  }

  if (!mesh_.is_distributed())
  {
    return;
  }

  //Build interior boundary
  boundary_ = new BoundaryMesh(mesh_, BoundaryMesh::interior);
  on_boundary_.init(mesh_, 0);
  for (VertexIterator vertex(*boundary_); !vertex.end(); ++vertex)
  {
    on_boundary_.set(boundary_->vertex_index(*vertex), true);
  }

  //
  uint gdim = mesh_.geometry().dim();
  MeshDistributedData& distdata = mesh_.distdata();

  _map<uint,std::vector<double> >::iterator owner_iterator=owner_tree.begin();
  _map<uint,std::vector<uint> >::iterator ghost_iterator=ghost_tree.begin();

  for (VertexIterator vertex(*boundary_); !vertex.end(); ++vertex)
  {
    Vertex on_mesh(mesh_, boundary_->vertex_index(*vertex));
    //Building owner tree:
    //The process number which owns the vertex is saved as key
    if (distdata[0].is_ghost(on_mesh.index()))
    {
      owner_iterator = owner_tree.find(distdata[0].get_owner(on_mesh.index()));
      if (owner_iterator != owner_tree.end())
      {
        (owner_iterator->second).push_back(
            double(distdata[0].get_global(on_mesh.index())));
      }
      else
      {
        std::vector<double> vertices_to_send;
        vertices_to_send.push_back(
            double(distdata[0].get_global(on_mesh.index())));
        owner_tree.insert(
            std::pair<uint, std::vector<double> >(
                distdata[0].get_owner(on_mesh.index()), vertices_to_send));
      }

      std::vector<double> vertex_info;
      double num_neigh = 0.0;
      double *sum = new double[gdim];
      for (uint j = 0; j < gdim; ++j)
      {
        sum[j] = 0.0;
      }
      std::vector<double> boundary_info;
      //building send_inner
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
      {
        if (on_mesh.index() == vn->index())
        {
          continue;
        }
        else
        {
          num_neigh += 1.0;
          // Compute center of mass
          const real* xn = vn->x();
          for (uint i = 0; i < gdim; ++i)
          {
            sum[i] += xn[i];
          }
        }
      }
      vertex_info.push_back(num_neigh);
      for (uint i = 0; i < gdim; ++i)
      {
        vertex_info.push_back(sum[i]);
      }
      send_inner.insert(
          std::pair<uint, std::vector<double> >(
              double(distdata[0].get_global(on_mesh.index())), vertex_info));
      delete[] sum;
    }

    //building recv_sum
    else
    {
      _set<uint> NeighboringProcessor = distdata[0].get_shared_adj(on_mesh.index());
      for (_set<uint>::iterator it = NeighboringProcessor.begin();
          it!= NeighboringProcessor.end();++it)
      {
        ghost_iterator=ghost_tree.find(*it);
        if(ghost_iterator!=ghost_tree.end())
        {
          (ghost_iterator->second).push_back(
              distdata[0].get_global(on_mesh.index()));
        }
        else
        {
          std::vector<uint> vertices_to_send;
          vertices_to_send.push_back(
              distdata[0].get_global(on_mesh.index()));
          ghost_tree.insert(
              std::pair<uint,std::vector<uint> >(*it, vertices_to_send));
        }
      }
      std::vector<double> vertex_info;
      std::vector<uint> vertex_check;
      double num_neigh = 0.0;
      double *sum = new double[gdim];

      for(uint j=0; j < gdim; ++j)
      {
        sum[j]=0.0;
      }
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
      {
        // Skip the vertex itself
        if (on_mesh.index() == vn->index())
        {
          continue;
        }
        num_neigh += 1.0;

        // Compute center of mass
        const real* xn = vn->x();
        for (uint i = 0; i < gdim; ++i)
        {
          sum[i] += xn[i];
        }
      }
      vertex_info.push_back(num_neigh);
      for (uint i = 0; i < gdim; ++i)
      {
        vertex_info.push_back(sum[i]);
      }
      recv_sum.insert(
          std::pair<uint,std::vector<double> >(
              distdata[0].get_global(on_mesh.index()), vertex_info));
      delete[] sum;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshSmoothData::sum_contribution(double*& recv_buff, int& mod,
                                      double& stopper, uint& src)
{

  if (!mesh_.is_distributed())
  {
    return;
  }

  int l = 0;
  _map<uint,std::vector<double> >::iterator receive_iterator=recv_sum.begin();

  while (recv_buff[l] != stopper)
  {
    receive_iterator = recv_sum.find(recv_buff[l]);
    if (receive_iterator != recv_sum.end())
    {
      for (uint j = 1; j <= (receive_iterator->second).size(); ++j)
      {
        (receive_iterator->second)[j - 1] += recv_buff[l + j];
      }
    }
    l += mod;
  }

}
//-----------------------------------------------------------------------------

}

