// Copyright (C) 2011 Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2011-06-30
// Last changed: 2011-06-30

#ifndef __MESH_SMOOTHDATA_H
#define __MESH_SMOOTHDATA_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <map>

namespace dolfin
{
  class MeshSmoothData
  {
    //    class Mesh;
    
  public:
    MeshSmoothData(Mesh& _mesh);
    ~MeshSmoothData();
    
    /*
      This function builds the map "owner_tree" which has the owner as
      key and attaches the number of neighbors, the sum in x direction,
      y direction and z direction. At the same time a map "recv_sum" is
      created with the same information but for the vertices owned by
      the core itself.x
    */
    void prepare_mesh();
    /*
      Summing up all the information which the process gets from the
      other process. In the map "ghost_tree" the source will be saved
      to able to send information back.
    */
    void sum_contribution(double*& recv_buff, 
			  int& mod, double& stopper, uint& src);
    
    _map<uint,std::vector<double> > owner_tree;
    _map<uint,std::vector<uint> > ghost_tree;
    _map<uint,std::vector<double> > send_inner;
    _map<uint,std::vector<double> > recv_sum;

    MeshFunction<bool> on_boundary;
    MeshFunction<bool> on_boundary_global;

    BoundaryMesh boundary;
    
  private:

    Mesh& mesh;    
  };  
}
#endif
