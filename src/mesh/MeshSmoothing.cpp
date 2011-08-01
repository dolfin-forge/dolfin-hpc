// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2011.
// Modified by Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller 2011.
//
// First added:  2008-07-16
// Last changed: 2011-06-30

#include <dolfin/common/constants.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshSmoothData.h>
#include <dolfin/mesh/MeshSmoothing.h>
#include <dolfin/parameter/parameters.h>
#include <map>
#include <vector>
#include <algorithm>


using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshSmoothing::smooth(Mesh& mesh)
{
  MeshSmoothData smooth_data(mesh);
  smooth_common(mesh, smooth_data);
}
//-----------------------------------------------------------------------------
void MeshSmoothing::smooth(Mesh& mesh, MeshSmoothData& smooth_data)
{
  smooth_common(mesh, smooth_data);
}
//-----------------------------------------------------------------------------
void MeshSmoothing::smooth_common(Mesh& mesh, MeshSmoothData& smooth_data)
{
   //starting MPI
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest, src;
  uint global_num_vertex=mesh.distdata().global_numVertices();
  double stopper= double(global_num_vertex+1);
    
  smooth_data.prepare_mesh();  
  // Create an local boundary mesh
  int boundary_number = (pe_size > 1 ? smooth_data.boundary().numVertices() : 0);
  const int d = mesh.geometry().dim();
  int module=d+2;//number of saved information, vertex index, number of neighbors, sum_x, sum_y, sum_z
  
  int max_un;
  int num_un = boundary_number * module;
#ifdef HAVE_MPI
    
  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INTEGER,MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  for (uint j = 1; j < pe_size; j++) 
  {
    double *recv_buff = new double[max_un];
    std::vector<double> send_buff;
    MPI_Status status;
    
    src = (rank -j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;
    _map<uint,std::vector<double> >::iterator sender_iterator=smooth_data.owner_tree.begin();
    sender_iterator=smooth_data.owner_tree.find(dest);
    
    if(sender_iterator!=smooth_data.owner_tree.end()){
      for (std::vector<double>::iterator iter_vector=sender_iterator->second.begin();
	   iter_vector != sender_iterator->second.end();iter_vector++)
      {
	send_buff.push_back(*iter_vector);//send key(as vertex index)
	std::vector<double> to_send_info;
	to_send_info=smooth_data.send_inner.find(*iter_vector)->second;
	
	for(uint i=0;i<to_send_info.size();i++){
	  send_buff.push_back(to_send_info[i]);
	}
      }
      
    }
    send_buff.push_back(stopper);
    MPI_Sendrecv(&send_buff[0], num_un, MPI_DOUBLE, dest, 1, recv_buff, max_un,
		 MPI_DOUBLE, src, 1,MPI::DOLFIN_COMM, &status);
    // sum inner contributions
    smooth_data.sum_contribution(recv_buff, module, stopper, src);

    delete[] recv_buff;
  }
  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
#endif
  
  // Iterate over all vertices
  std::vector<real> xx(d);
  _map<uint,std::vector<double> >::iterator boundary_iterator=smooth_data.recv_sum.begin();
  _map<uint,std::vector<double> >::iterator receive_iterator=smooth_data.recv_sum.begin();

  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    
    // Get coordinates of vertex
    real* x = v->x();
    const Point p = v->point();
    
    // Skip vertices on the boundary 
    double num_neighbors = 0.0;
    if(smooth_data.on_boundary_global()(*v))
      continue;
    else if(smooth_data.on_boundary()(*v) && pe_size > 1){
      receive_iterator=smooth_data.recv_sum.find(mesh.distdata().get_global(v->index(),0));
      if(receive_iterator!=smooth_data.recv_sum.end()){
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
    else
    {
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
    
    if(dolfin_get("Mesh smoothing restricted by rmin")) 
    {
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
    else
      for (uint i = 0; i < d; i++)
	x[i] =xx[i];
  }
  
#ifdef HAVE_MPI
  //sending information back
  for (uint j = 1; j < pe_size; j++) 
  {
    double *recv_buff = new double[max_un];
    std::vector<double> send_buff;
    MPI_Status status;
    
    src = (rank -j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;
    _map<uint,std::vector<uint> >::iterator sender_iterator= smooth_data.ghost_tree.begin(); //key=processor rank, 
    if(smooth_data.ghost_tree.size()!=0){
      sender_iterator=smooth_data.ghost_tree.find(dest);
      
      if(sender_iterator!=smooth_data.ghost_tree.end()){
	for (std::vector<uint>::iterator iter_vector=sender_iterator->second.begin();
	     iter_vector != sender_iterator->second.end();iter_vector++)
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
#endif
}
//-----------------------------------------------------------------------------
