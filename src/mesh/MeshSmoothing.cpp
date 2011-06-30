// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2010.
// Modified by Jeanentte Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller 2011.
// First added:  2008-07-16
// Last changed: 2011-04-23

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include "MeshSmoothing.h"
#include <mpi.h>
#include <map>
#include <vector>
#include <algorithm>
#include <dolfin/main/MPI.h>


using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshSmoothing::prepare_mesh(std::map<uint,std::vector<double> >& owner_tree, 
				 std::map<uint,std::vector<uint> >& ghost_tree,
				 std::map<uint,std::vector<double> >& send_inner,
				 std::map<uint,std::vector<double> >& recv_sum,
				 BoundaryMesh& boundary, 
				 Mesh& mesh, MeshFunction<uint>*& vertex_map,
				 MeshFunction<bool>& on_boundary,int d)
{
  std::map<uint,std::vector<double> >::iterator owner_iterator=owner_tree.begin();
  std::map<uint,std::vector<uint> >::iterator ghost_iterator=ghost_tree.begin();

  for (VertexIterator vertex(boundary); !vertex.end(); ++vertex){
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
      double *sum = new double[d];
      for(int j=0;j<d;j++)
	sum[j]=0.0;
      std::vector<double> boundary_info;
      //building send_inner
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
	{
	  if (on_mesh.index() == vn->index())
	    continue;
	  //else if(!on_boundary.get(vn->index())){
	  else{
	    num_neigh += 1.0;
	    // Compute center of mass
	    const real* xn = vn->x();
	    for (int i = 0; i < d; i++)
	      sum[i] += xn[i];
	  }
	}
      vertex_info.push_back(num_neigh); 
      for (int i = 0; i < d; i++){
	vertex_info.push_back(sum[i]);
      }
      send_inner.insert(std::pair<uint,std::vector<double> >(double(mesh.distdata().get_global(on_mesh.index(),0)), vertex_info));
      delete[] sum;
    }
    
    //building recv_sum
    else{
      _set<uint> NeighboringProcessor = mesh.distdata().get_shared_adj(on_mesh.index(),0);
      for (_set<uint>::iterator it =  NeighboringProcessor.begin(); it!= NeighboringProcessor.end();it++){
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
      double *sum = new double[d];
      
      for(int j=0;j<d;j++)
	sum[j]=0.0;
      for (VertexIterator vn(on_mesh); !vn.end(); ++vn)
	{
	  // Skip the vertex itself
	  if (on_mesh.index() == vn->index())
	    continue;
	  num_neigh += 1.0;
	  
	  // Compute center of mass
	  const real* xn = vn->x();
	  for (int i = 0; i < d; i++)
	    sum[i] += xn[i];
	}
      vertex_info.push_back(num_neigh); 
      for (int i = 0; i < d; i++){
	vertex_info.push_back(sum[i]);
      }
      recv_sum.insert(std::pair<uint,std::vector<double> >(mesh.distdata().get_global(on_mesh.index(),0), vertex_info));
      delete[] sum;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshSmoothing::sum_contribution(std::map<uint,std::vector<double> >& recv_sum,
				     double*& recv_buff, 
				     int& mod, double& stopper, uint& src)
{
  
  int l=0;
  std::map<uint,std::vector<double> >::iterator receive_iterator=recv_sum.begin();
  //std::map<uint,std::vector<uint> >::iterator ghost_iterator=ghost_tree.begin();
  
  while(recv_buff[l]!=stopper){
    receive_iterator=recv_sum.find(recv_buff[l]);
    if(receive_iterator!=recv_sum.end())
      for (uint j=1;j<=(receive_iterator->second).size();j++){
	(receive_iterator->second)[j-1]+=recv_buff[l+j];
      }
    /*
    ghost_iterator=ghost_tree.find(src);
    if(ghost_iterator!=ghost_tree.end()){
      (ghost_iterator->second).push_back(recv_buff[l]);
    }
    else{
      std::vector<uint> ghost_info;
      ghost_info.push_back(recv_buff[l]);
      ghost_tree.insert(std::pair<uint,std::vector<uint> >(src, ghost_info));
    }
    */
    l+=mod;
  }
  
}
//-----------------------------------------------------------------------------
void MeshSmoothing::smooth(Mesh& mesh)
{
   //starting MPI
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest, src;
  uint global_num_vertex=mesh.distdata().global_numVertices();
  double stopper= double(global_num_vertex+1);
  
  std::map<uint,std::vector<double> > owner_tree; //saves owner as key
  //to trace sender
  std::map<uint,std::vector<uint> > ghost_tree; //key=processor rank, 
  //send maps
  std::map<uint,std::vector<double> > send_inner; //saves information about inner incident nodes, key=vertex index
  //receive maps
  std::map<uint,std::vector<double> > recv_sum; //sum up x, y, z, and number of neighbors, key=vertex index
  
  //mapping for different boundaries
  //interior
  BoundaryMesh boundary;
  boundary.init_interior(mesh);
  std::cout << "boudary size:  " <<boundary.numVertices()<<std::endl;
  message("Rank: %d has %d ghosted and %d shared vertices", dolfin::MPI::processNumber(), mesh.distdata().num_ghost(0), mesh.distdata().num_shared(0));
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");
  dolfin_assert(vertex_map);
  MeshFunction<bool> on_boundary(mesh,0);
  on_boundary=false;
  for (VertexIterator vertex(boundary); !vertex.end(); ++vertex){
     on_boundary.set((*vertex_map)(*vertex), true);
  }
  
  //global boundary
  BoundaryMesh boundary_global;
  boundary_global.init(mesh);
  MeshFunction<uint>* vertex_map_global = boundary_global.data().meshFunction("vertex map");
  dolfin_assert(vertex_map_global);
  MeshFunction<bool> on_boundary_global(mesh, 0);
  on_boundary_global=false;
  if(boundary_global.numVertices()!=0){
    for (VertexIterator vertex(boundary_global); !vertex.end(); ++vertex){
      on_boundary_global.set((*vertex_map_global)(*vertex), true);
    }
  }
  
  // Create an local boundary mesh
  int boundary_number=boundary.numVertices();
  const int d = mesh.geometry().dim();
  int module=d+2;//number of saved information, vertex index, number of neighbors, sum_x, sum_y, sum_z
  MeshSmoothing::prepare_mesh(owner_tree, ghost_tree,send_inner, recv_sum, boundary, mesh, vertex_map, on_boundary,d);
  int max_un;
  int num_un = boundary_number * module;
 
  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INTEGER,MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  
  for (uint j = 1; j < pe_size; j++) 
    {
      double *recv_buff = new double[max_un];
      Array<double> send_buff;
      MPI_Status status;
            
      src = (rank -j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;
      std::map<uint,std::vector<double> >::iterator sender_iterator=owner_tree.begin();
      sender_iterator=owner_tree.find(dest);
      
      if(sender_iterator!=owner_tree.end()){
	for (std::vector<double>::iterator iter_vector=sender_iterator->second.begin();iter_vector != sender_iterator->second.end();iter_vector++)
	  {
	    send_buff.push_back(*iter_vector);//send key(as vertex index)
	    std::vector<double> to_send_info;
	    to_send_info=send_inner.find(*iter_vector)->second;
	    
	    for(uint i=0;i<to_send_info.size();i++){
	      send_buff.push_back(to_send_info[i]);
	    }
	  }
      
      }
      send_buff.push_back(stopper);
      MPI_Sendrecv(&send_buff[0], num_un, MPI_DOUBLE, dest, 1, recv_buff, max_un, MPI_DOUBLE, src, 1,MPI::DOLFIN_COMM, &status);
      // sum inner contributions
      MeshSmoothing::sum_contribution(recv_sum,recv_buff, module, stopper, src);
      delete[] recv_buff;
    }

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);

  // Iterate over all vertices

  Array<real> xx(d);
  std::map<uint,std::vector<double> >::iterator boundary_iterator=recv_sum.begin();
  std::map<uint,std::vector<double> >::iterator receive_iterator=recv_sum.begin();

  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    
    // Get coordinates of vertex
    real* x = v->x();
    const Point p = v->point();
    
    // Skip vertices on the boundary 
    double num_neighbors = 0.0;
    if(on_boundary_global(*v))
      continue;
    else if(on_boundary(*v)){
      receive_iterator=recv_sum.find(mesh.distdata().get_global(v->index(),0));
      if(receive_iterator!=recv_sum.end()){
	for (int i = 0; i < d; i++) xx[i] = 0.0;
	num_neighbors = (receive_iterator->second)[0];
	for (int i = 0; i < d; i++)
	  {
	    xx[i] += (receive_iterator->second)[i+1];
	  }
	
	for (int i = 0; i < d; i++)
	  xx[i] /= static_cast<real>(num_neighbors); 
      }
    }

    else{
      // Compute center of mass of neighboring vertices
      for (int i = 0; i < d; i++) 
	xx[i] = 0.0;
      num_neighbors = 0.0;
      
      for (VertexIterator vn(*v); !vn.end(); ++vn)
	{
	  // Skip the vertex itself
	  if (v->index() == vn->index())
	    continue;
	  
	  num_neighbors += 1.0;
	  
	  // Compute center of mass
	  const real* xn = vn->x();
	  for (int i = 0; i < d; i++)
	    xx[i] += xn[i];
	}
      
      for (int i = 0; i < d; i++)
	xx[i] /= static_cast<real>(num_neighbors);
      
    }

    // Compute closest distance to boundary of star
    real rmin = 0.0;
    for (CellIterator c(*v); !c.end(); ++c)
    {
      // Get local number of vertex relative to facet
      const uint local_vertex = c->index(*v);

      // Get normal of corresponding facet
      Point n = c->normal(local_vertex);
      
      // Get first vertex in facet
      Facet f(mesh, c->entities(mesh.topology().dim() - 1)[local_vertex]);
      VertexIterator fv(f);

      // Compute length of projection of v - fv onto normal
      const real r = std::abs(n.dot(p - fv->point()));
      if (rmin == 0.0)
        rmin = r;
      else
        rmin = std::min(rmin, r);
    }

    // Move vertex at most a distance rmin / 2
    real r = 0.0;
    for (int i = 0; i < d; i++)
    {
      const real dx = xx[i] - x[i];
      r += dx*dx;
    }
    r = std::sqrt(r);
    if (r < DOLFIN_EPS)
      continue;
    rmin = std::min(0.5*rmin, r);
    for (int i = 0; i < d; i++)
      x[i] += rmin*(xx[i] - x[i])/r;
  }

  //sending information back
  for (uint j = 1; j < pe_size; j++) 
    {
      double *recv_buff = new double[max_un];
      Array<double> send_buff;
      MPI_Status status;
      
      src = (rank -j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;
      std::map<uint,std::vector<uint> >::iterator sender_iterator= ghost_tree.begin(); //key=processor rank, 
      if(ghost_tree.size()!=0){
	sender_iterator=ghost_tree.find(dest);
	
	if(sender_iterator!=ghost_tree.end()){
	  for (std::vector<uint>::iterator iter_vector=sender_iterator->second.begin();iter_vector != sender_iterator->second.end();iter_vector++)
	    {
	      send_buff.push_back(*iter_vector);// global vertex index
	      std::vector<double> to_send_info;
	      Vertex on_mesh(mesh, mesh.distdata().get_local((*iter_vector),0));
	      for(int i=0;i<d;i++){
		send_buff.push_back(on_mesh.x(i));//send x, y, z
	      }
	    }
	}
      }
      send_buff.push_back(stopper);
      MPI_Sendrecv(&send_buff[0], num_un, MPI_DOUBLE, dest, 1, recv_buff, max_un, MPI_DOUBLE, src, 1,MPI::DOLFIN_COMM, &status);
      
      int l=0;
      while(recv_buff[l]!=stopper){
	Vertex on_mesh(mesh, mesh.distdata().get_local(recv_buff[l],0));
	// Get coordinates of vertex
	real* x = on_mesh.x();
	for(int j=0; j<d;j++){
	  x[j]=recv_buff[l+j+1];
	}
	l+=d+1;
      }
      delete[] recv_buff;
    }
  
  //MPI_Finalize();
}
//-----------------------------------------------------------------------------
