// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_DEPENDENT_H
#define __DOLFIN_MESH_DEPENDENT_H

#include <dolfin/mesh/Mesh.h>

#include <string>

namespace dolfin
{

class MeshDependent
{

public:

  ///
  explicit MeshDependent(Mesh& mesh);

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

  ///
  MeshDependent(MeshDependent const& other);

  Mesh * const mesh_;
  int topology_token_;
  int geometry_token_;

};

//---------------------------------------------------------------------------
inline Mesh& MeshDependent::mesh() const
{
  return *mesh_;
}

//---------------------------------------------------------------------------
inline bool MeshDependent::invalid_mesh_topology() const
{
  return topology_token_ != mesh_->topology().token();
}

//---------------------------------------------------------------------------
inline bool MeshDependent::invalid_mesh_geometry() const
{
  return geometry_token_ != mesh_->geometry().token();
}

//---------------------------------------------------------------------------
inline bool MeshDependent::invalid_mesh() const
{
  return invalid_mesh_topology() || invalid_mesh_geometry();
}

//---------------------------------------------------------------------------
inline void MeshDependent::update_mesh_dependency()
{
  topology_token_ = mesh_->topology().token();
  geometry_token_ = mesh_->geometry().token();
}

}

#endif /* __DOLFIN_MESH_DEPENDENT_H */
