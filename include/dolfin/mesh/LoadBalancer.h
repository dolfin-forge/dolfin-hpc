// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LOAD_BALANCER_H
#define __DOLFIN_LOAD_BALANCER_H

#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Edge;
class Vertex;

class LoadBalancer
{

public:
  using VertexFunction      = MeshValues< real, Vertex >;
  using CellFunction        = MeshValues< size_t, Cell >;
  using VertexFunctionPPair = std::pair< VertexFunction *, VertexFunction * >;
  using CellFunctionPPair   = std::pair< CellFunction *, CellFunction * >;

  /// Balancing types
  enum Type
  {
    Default,     //< default
    LEPP,        //< longest edge propagation paths (for RivaraRefinement)
    EdgeCollapse //< optimized for edge collapse
  };

  /// Balance mesh according to predefined weight function
  static void balance( Mesh & mesh, MeshValues< size_t, Cell > & weight );

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh
  static void balance( Mesh &                     mesh,
                       MeshValues< bool, Cell > & cell_marker,
                       Type                       type = Default );

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance( Mesh &                     mesh,
                       MeshValues< bool, Cell > & cell_marker,
                       real                       tf,
                       real                       tb,
                       real                       ts,
                       Type                       type = Default );

  /// Balanace mesh according to predefined weight function and preserve mesh
  /// functions
  static void balance( Mesh &                               mesh,
                       MeshValues< size_t, Cell > &         weight,
                       std::vector< VertexFunctionPPair > & vertex_functions );

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh
  /// functions
  static void balance( Mesh &                               mesh,
                       MeshValues< bool, Cell > &           cell_marker,
                       std::vector< VertexFunctionPPair > & vertex_functions,
                       Type                                 type = Default );

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance( Mesh &                               mesh,
                       MeshValues< bool, Cell > &           cell_marker,
                       std::vector< VertexFunctionPPair > & vertex_functions,
                       real                                 tf,
                       real                                 tb,
                       real                                 ts,
                       Type                                 type = Default );

  /// Balanace mesh according to predefined weight function and preserve mesh
  /// functions (vertex and cell)
  static void balance( Mesh &                               mesh,
                       MeshValues< size_t, Cell > &         weight,
                       std::vector< CellFunctionPPair > &   cell_functions,
                       std::vector< VertexFunctionPPair > & vertex_functions );

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh
  /// functions(vertex and cell)
  static void balance( Mesh &                               mesh,
                       MeshValues< bool, Cell > &           cell_marker,
                       std::vector< CellFunctionPPair > &   cell_functions,
                       std::vector< VertexFunctionPPair > & vertex_functions,
                       Type                                 type = Default );

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions (vertex and
  /// cell) tb    Time to perform one flop in seconds tb    Time to transmit one
  /// byte in seconds ts    Startup time in seconds (Interconnect latency)
  static void balance( Mesh &                               mesh,
                       MeshValues< bool, Cell > &           cell_marker,
                       std::vector< CellFunctionPPair > &   cell_functions,
                       std::vector< VertexFunctionPPair > & vertex_functions,
                       real                                 tf,
                       real                                 tb,
                       real                                 ts,
                       Type                                 type = Default );

  /// Balanace mesh according to predefined weight function and preserve mesh
  /// functions (cell)
  static void balance( Mesh &                             mesh,
                       MeshValues< size_t, Cell > &       weight,
                       std::vector< CellFunctionPPair > & cell_functions );

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh
  /// functions (cell)
  static void balance( Mesh &                             mesh,
                       MeshValues< bool, Cell > &         cell_marker,
                       std::vector< CellFunctionPPair > & cell_functions,
                       Type                               type = Default );

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions (cell)
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance( Mesh &                             mesh,
                       MeshValues< bool, Cell > &         cell_marker,
                       std::vector< CellFunctionPPair > & cell_functions,
                       real                               tf,
                       real                               tb,
                       real                               ts,
                       Type                               type = Default );

private:
  static void weight_function( Mesh &                       mesh,
                               MeshValues< bool, Cell > &   cell_marker,
                               MeshValues< size_t, Cell > & weight,
                               size_t *                     w_sum,
                               Type                         type );

  static void weight_lepp( Mesh &                       mesh,
                           Cell &                       c,
                           Edge &                       ce,
                           MeshValues< size_t, Cell > & weight,
                           size_t                       depth );

  static void process_reassignment( MeshFunction< size_t > & partitions,
                                    size_t *                 max_sendrecv );

  static auto computational_gain( Mesh &                       mesh,
                                  MeshValues< size_t, Cell > & weight,
                                  MeshValues< size_t, Cell > & partitions,
                                  size_t                       max_sendrecv,
                                  real                         tf,
                                  real                         tb,
                                  real                         ts ) -> bool;

  static void
    radixsort_matrix( size_t * res, size_t * Matrix, size_t m, bool desc );

  static void pradixsort_matrix( size_t * res, size_t * Matrix, size_t m );

public:
  // Very imperfect but good enough for now
  static auto partitions( Mesh & mesh ) -> MeshValues< size_t, Cell > &;

  static auto clear( Mesh & mesh ) -> bool;

private:
  static _ordered_map< Mesh *, MeshValues< size_t, Cell > * > s_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_LOAD_BALANCER_H */
