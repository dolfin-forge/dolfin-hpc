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
  virtual auto mesh() const -> Mesh&;

  ///
  auto invalid_mesh_topology() const -> bool;

  ///
  auto invalid_mesh_geometry() const -> bool;

  ///
  auto invalid_mesh() const -> bool;

  ///
  void update_mesh_dependency();

protected:

  ///
  virtual ~MeshDependent() = default;

private:

  ///
  MeshDependent(MeshDependent const& other);

  Mesh * const mesh_;
  int topology_token_;
  int geometry_token_;

};

//---------------------------------------------------------------------------
inline auto MeshDependent::mesh() const -> Mesh&
{
  return *mesh_;
}

//---------------------------------------------------------------------------
inline auto MeshDependent::invalid_mesh_topology() const -> bool
{
  return topology_token_ != mesh_->topology().token();
}

//---------------------------------------------------------------------------
inline auto MeshDependent::invalid_mesh_geometry() const -> bool
{
  return geometry_token_ != mesh_->geometry().token();
}

//---------------------------------------------------------------------------
inline auto MeshDependent::invalid_mesh() const -> bool
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
