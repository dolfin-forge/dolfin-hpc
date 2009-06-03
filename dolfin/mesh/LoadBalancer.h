// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//
//
// First added:  2008-03-03
// Last changed: 2009-05-24

#ifndef __LOAD_BALANCER_H
#define __LOAD_BALANCER_H

#include "MeshFunction.h"

namespace dolfin
{
  class LoadBalancer
  {
  public:

    enum Type { Default, LEPP};

    // Balanace mesh according to predefined weight function
    static void balance(Mesh& mesh, MeshFunction<uint>& weight);

    // Balance mesh according to marked cells, 
    // new_cell_marker marks cells in cell_marker for new mesh
    static void balance(Mesh& mesh, MeshFunction<bool>& cell_marker, 
			Type type = Default);

    // Balance mesh according to marked cells, tune loadbalancer
    // with machine specific parameters
    // tb    Time to perform one flop in seconds
    // tb    Time to transmit one byte in seconds
    // ts    Startup time in seconds (Interconnect latency)
    static void balance(Mesh& mesh, MeshFunction<bool>& cell_marker,
			real tf, real tb, real ts, Type type = Default);

  private:

    static void weight_function(Mesh& mesh, 
				MeshFunction<bool>& cell_marker,
				MeshFunction<uint>& weight,
				uint* w_sum, Type type);

    static void weight_lepp(Mesh& mesh, Cell& c, Edge& ce,
			    MeshFunction<uint>& weight, uint depth);

    static void process_reassignment(MeshFunction<uint>& partitions,
				     uint* max_sendrecv);

    static bool computational_gain(Mesh& mesh, 
				   MeshFunction<uint>& weight,
				   MeshFunction<uint>& partitions,
				   uint max_sendrecv,
				   real tf, real tb, real ts);

    static void radixsort_matrix(uint* res, uint* Matrix, uint m, bool desc);

    static void pradixsort_matrix(uint* res, uint* Matrix, uint m);
  };
}

#endif


