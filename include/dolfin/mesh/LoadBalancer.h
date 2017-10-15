// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Balthasar Reuter, 2013.
//
// First added:  2008-03-03
// Last changed: 2013-05-31

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

  typedef MeshValues<real, Vertex> VertexFunction;
  typedef MeshValues<uint, Cell> CellFunction;
  typedef std::pair<VertexFunction *, VertexFunction *> VertexFunctionPPair;
  typedef std::pair<CellFunction *, CellFunction *> CellFunctionPPair;

  /// Balancing types
  enum Type
  {
    Default,    //< default
    LEPP,       //< longest edge propagation paths (for RivaraRefinement)
    EdgeCollapse //< optimized for edge collapse
  };

  /// Balance mesh according to predefined weight function
  static void balance(Mesh& mesh, MeshValues<uint, Cell>& weight);

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Type type = Default);

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker, real tf,
                      real tb, real ts, Type type = Default);

  /// Balanace mesh according to predefined weight function and preserve mesh functions
  static void balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                      Array<VertexFunctionPPair>& vertex_functions);

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh functions
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<VertexFunctionPPair>& vertex_functions,
                      Type type = Default);

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<VertexFunctionPPair>& vertex_functions, real tf,
                      real tb, real ts, Type type = Default);

  /// Balanace mesh according to predefined weight function and preserve mesh functions (vertex and cell)
  static void balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                      Array<CellFunctionPPair>& cell_functions,
                      Array<VertexFunctionPPair>& vertex_functions);

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh functions(vertex and cell)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<CellFunctionPPair>& cell_functions,
                      Array<VertexFunctionPPair>& vertex_functions,
                      Type type = Default);

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions (vertex and cell)
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<CellFunctionPPair>& cell_functions,
                      Array<VertexFunctionPPair>& vertex_functions, real tf,
                      real tb, real ts, Type type = Default);

  /// Balanace mesh according to predefined weight function and preserve mesh functions (cell)
  static void balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                      Array<CellFunctionPPair>& cell_functions);

  /// Balance mesh according to marked cells,
  /// new_cell_marker marks cells in cell_marker for new mesh and  preserve mesh functions (cell)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<CellFunctionPPair>& cell_functions,
                      Type type = Default);

  /// Balance mesh according to marked cells, tune loadbalancer
  /// with machine specific parameters and  preserve mesh functions (cell)
  /// tb    Time to perform one flop in seconds
  /// tb    Time to transmit one byte in seconds
  /// ts    Startup time in seconds (Interconnect latency)
  static void balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                      Array<CellFunctionPPair>& cell_functions, real tf,
                      real tb, real ts, Type type = Default);

private:

  static void weight_function(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                              MeshValues<uint, Cell>& weight, uint* w_sum,
                              Type type);

  static void weight_lepp(Mesh& mesh, Cell& c, Edge& ce,
                          MeshValues<uint, Cell>& weight, uint depth);

  static void process_reassignment(MeshFunction<uint>& partitions,
                                   uint* max_sendrecv);

  static bool computational_gain(Mesh& mesh,
                                 MeshValues<uint, Cell>& weight,
                                 MeshValues<uint, Cell>& partitions,
                                 uint max_sendrecv, real tf, real tb, real ts);

  static void radixsort_matrix(uint* res, uint* Matrix, uint m, bool desc);

  static void pradixsort_matrix(uint* res, uint* Matrix, uint m);

};

} /* namespace dolfin */

#endif /* __DOLFIN_LOAD_BALANCER_H */
