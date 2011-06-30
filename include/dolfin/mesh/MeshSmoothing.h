// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Jeanentte Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller 2011.
// First added:  2008-07-16
// Last changed: 2011-06-30

#ifndef __MESH_SMOOTHING_H
#define __MESH_SMOOTHING_H

#include <dolfin/mesh/BoundaryMesh.h>
#include <map>

namespace dolfin
{
  
  class Mesh;

  /// This class implements mesh smoothing. The coordinates of
  /// internal vertices are updated by local averaging.

  class MeshSmoothing
  {

  public:

    static void smooth(Mesh& mesh);

    /*
    This function builds the map "owner_tree" which has the owner as
    key and attaches the number of neighbors, the sum in x direction,
    y direction and z direction. At the same time a map "recv_sum" is
    created with the same information but for the vertices owned by
    the core itself.x
    */
    static void prepare_mesh(std::map<uint,std::vector<double> >& owner_tree,
			     std::map<uint,std::vector<uint> >& ghost_tree,
			     std::map<uint,std::vector<double> >& send_inner,
			     std::map<uint,std::vector<double> >& recv_sum,
			     BoundaryMesh& boundary, Mesh& mesh,
			     MeshFunction<uint>*& vertex_map,
			     MeshFunction<bool>& on_boundary,int d);
    /*
      Summing up all the information which the process gets from the
      other process. In the map "ghost_tree" the source will be saved
      to able to send information back.
    */
    static void sum_contribution(std::map<uint, std::vector<double> >& recv_sum, 
				 double*& recv_buff, 
				 int& mod, double& stopper, uint& src);
    
  };

}

#endif
