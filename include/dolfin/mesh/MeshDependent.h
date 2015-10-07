// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-04
// Last changed: 2014-03-04

#ifndef __DOLFIN_MESH_DEPENDENT_H
#define __DOLFIN_MESH_DEPENDENT_H

#include <string>

namespace dolfin
{

class Mesh;

class MeshDependent
{
public:

  ///
  MeshDependent(Mesh& mesh);

  ///
  virtual Mesh& mesh() const;

  ///
  bool invalid_mesh_topology() const;

  ///
  bool invalid_mesh_geometry() const;

  ///
  bool invalid_mesh() const;

  ///
  void update_mesh_dependency();

protected:

  ///
  virtual ~MeshDependent();

private:

  Mesh& mesh_;
  int topology_token_;
  int geometry_token_;

};

}

#endif /* __DOLFIN_MESH_DEPENDENT_H */
