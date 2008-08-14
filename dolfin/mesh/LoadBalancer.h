// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//
//
// First added:  2008-03-03
// Last changed: 2008-07-17

#ifndef __LOAD_BALANCER_H
#define __LOAD_BALANCER_H

#include "MeshFunction.h"
#include <dolfin/main/MPI.h>

namespace dolfin
{
  class LoadBalancer
  {
  public:

    // Balance mesh according to marked cells, 
    // new_cell_marker marks cells in cell_marker for new mesh
    static void balance(Mesh& mesh, MeshFunction<bool>& cell_marker);

    // Balance mesh according to marked cells, tune loadbalancer
    // with machine specific parameters
    // tb    Time to perform one flop in seconds
    // tb    Time to transmit one byte in seconds
    // ts    Startup time in seconds (Interconnect latency)
    static void balance(Mesh& mesh, MeshFunction<bool>& cell_marker,
			real tf, real tb, real ts);

  private:

    static void weight_function(Mesh& mesh, 
				MeshFunction<bool>& cell_marker,
				MeshFunction<uint>& weight,
				uint* w_sum);

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


