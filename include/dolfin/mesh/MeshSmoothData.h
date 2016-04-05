// Copyright (C) 2011 Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2011-06-30
// Last changed: 2011-06-30

#ifndef __DOLFIN_MESH_SMOOTHDATA_H
#define __DOLFIN_MESH_SMOOTHDATA_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <map>

namespace dolfin
{

class MeshSmoothData
{

public:

  ///
  MeshSmoothData(Mesh& mesh);

  ///
  ~MeshSmoothData();

  /**
   This function builds the map "owner_tree" which has the owner as
   key and attaches the number of neighbors, the sum in x direction,
   y direction and z direction. At the same time a map "recv_sum" is
   created with the same information but for the vertices owned by
   the core itself.x
   */
  void prepare_mesh();

  /**
   Summing up all the information which the process gets from the
   other process. In the map "ghost_tree" the source will be saved
   to able to send information back.
   */
  void sum_contribution(real*& recv_buff, int& mod, real& stopper, uint& src);

  ///
  BoundaryMesh& boundary();

  ///
  MeshFunction<bool>& on_boundary();

  ///
  MeshFunction<bool>& on_boundary_global();

  ///
  _map<uint,std::vector<real> > owner_tree;
  _map<uint,std::vector<uint> > ghost_tree;
  _map<uint,std::vector<real> > send_inner;
  _map<uint,std::vector<real> > recv_sum;

private:

  Mesh& mesh_;
  BoundaryMesh * boundary_;

  MeshFunction<bool> on_boundary_;
  MeshFunction<bool> on_boundary_global_;
};

//--- INLINES -----------------------------------------------------------------

inline BoundaryMesh& MeshSmoothData::boundary()
{
  return *boundary_;
}

//-----------------------------------------------------------------------------
inline MeshFunction<bool>& MeshSmoothData::on_boundary()
{
  return on_boundary_;
}

//-----------------------------------------------------------------------------
inline MeshFunction<bool>& MeshSmoothData::on_boundary_global()
{
  return on_boundary_global_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SMOOTHDATA_H */
