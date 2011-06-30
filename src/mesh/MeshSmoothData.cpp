// Copyright (C) 2011 Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2011-06-30
// Last changed: 2011-06-30

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/MeshSmoothData.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshSmoothData::MeshSmoothData(Mesh& _mesh) : mesh(_mesh)
{


}
//-----------------------------------------------------------------------------
MeshSmoothData::~MeshSmoothData()
{
}
//-----------------------------------------------------------------------------
void MeshSmoothData::prepare_mesh()
{
  //mapping for different boundaries
  on_boundary.init(mesh, 0);
  on_boundary = false;
  
  on_boundary_global.init(mesh, 0);
  on_boundary_global = false;

  //interior
  _boundary.init_interior(mesh);
  MeshFunction<uint>* vertex_map = _boundary.data().meshFunction("vertex map");
  dolfin_assert(vertex_map);

  for (VertexIterator vertex(_boundary); !vertex.end(); ++vertex)
     on_boundary.set((*vertex_map)(*vertex), true);
  
  //global boundary
  BoundaryMesh boundary_global;
  boundary_global.init(mesh);
  MeshFunction<uint>* vertex_map_global = boundary_global.data().meshFunction("vertex map");
  dolfin_assert(vertex_map_global);
  on_boundary_global=false;
  if(boundary_global.numVertices()!=0){
    for (VertexIterator vertex(boundary_global); !vertex.end(); ++vertex){
      on_boundary_global.set((*vertex_map_global)(*vertex), true);
    }
  }
  uint gdim = mesh.geometry().dim();
 
  _map<uint,std::vector<double> >::iterator owner_iterator=owner_tree.begin();
  _map<uint,std::vector<uint> >::iterator ghost_iterator=ghost_tree.begin();

  for (VertexIterator vertex(_boundary); !vertex.end(); ++vertex){
    Vertex on_mesh(mesh,  vertex_map->get(*vertex));
    //building owner tree: the number of the CPU which owns the vertex is saved as key 
    if(mesh.distdata().is_ghost(on_mesh.index(), 0)){
      owner_iterator=owner_tree.find(mesh.distdata().get_owner(on_mesh.index(),0));
      if(owner_iterator!=owner_tree.end())
	(owner_iterator->second).push_back(double(mesh.distdata().get_global(on_mesh.index(),0)));
      else
	{
	  std::vector<double> vertices_to_send;
	  vertices_to_send.push_back(double(mesh.distdata().get_global(on_mesh.index(),0)));
	  owner_tree.insert(std::pair<uint,std::vector<double> >(mesh.distdata().get_owner(on_mesh.index(),0), vertices_to_send));
	}
      
      std::vector<double> vertex_info;
      double num_neigh = 0.0;
      double *sum = new double[gdim];
      for(int j=0; j < gdim; j++)
	sum[j]=0.0;
      std::vector<double> boundary_info;
      //building send_inner
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
      {
	if (on_mesh.index() == vn->index())
	  continue;
	else{
	  num_neigh += 1.0;
	  // Compute center of mass
	  const real* xn = vn->x();
	  for (int i = 0; i < gdim; i++)
	    sum[i] += xn[i];
	  }
      }
      vertex_info.push_back(num_neigh); 
      for (int i = 0; i < gdim; i++)
      {
	vertex_info.push_back(sum[i]);
      }
      send_inner.insert(std::pair<uint,std::vector<double> >(double(mesh.distdata().get_global(on_mesh.index(),0)), vertex_info));
      delete[] sum;
    }
    
    //building recv_sum
    else{
      _set<uint> NeighboringProcessor = mesh.distdata().get_shared_adj(on_mesh.index(),0);
      for (_set<uint>::iterator it =  NeighboringProcessor.begin(); it!= NeighboringProcessor.end();it++)
      {
	ghost_iterator=ghost_tree.find(*it);
	if(ghost_iterator!=ghost_tree.end())
	  (ghost_iterator->second).push_back(mesh.distdata().get_global(on_mesh.index(),0));
	else
	{
	  std::vector<uint> vertices_to_send;
	  vertices_to_send.push_back(mesh.distdata().get_global(on_mesh.index(),0));
	  ghost_tree.insert(std::pair<uint,std::vector<uint> >(*it, vertices_to_send));
	}
      }
      std::vector<double> vertex_info;
      std::vector<uint> vertex_check;
      double num_neigh = 0.0;
      double *sum = new double[gdim];
      
      for(int j=0; j < gdim; j++)
	sum[j]=0.0;
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
      {
	// Skip the vertex itself
	if (on_mesh.index() == vn->index())
	  continue;
	num_neigh += 1.0;
	
	// Compute center of mass
	const real* xn = vn->x();
	for (int i = 0; i < gdim; i++)
	  sum[i] += xn[i];
      }
      vertex_info.push_back(num_neigh); 
      for (int i = 0; i < gdim; i++)
	vertex_info.push_back(sum[i]);
      recv_sum.insert(std::pair<uint,std::vector<double> >(mesh.distdata().get_global(on_mesh.index(),0), vertex_info));
      delete[] sum;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshSmoothData::sum_contribution(double*& recv_buff, int& mod,
				      double& stopper, uint& src)
{
  
  int l=0;
  _map<uint,std::vector<double> >::iterator receive_iterator=recv_sum.begin();
  
  while(recv_buff[l]!=stopper){
    receive_iterator=recv_sum.find(recv_buff[l]);
    if(receive_iterator!=recv_sum.end())
      for (uint j=1;j<=(receive_iterator->second).size();j++){
	(receive_iterator->second)[j-1]+=recv_buff[l+j];
      }
    l+=mod;
  }
  
}
//-----------------------------------------------------------------------------
