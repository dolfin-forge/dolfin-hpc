// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#ifndef __PTESTS_MESH_H_
#define __PTESTS_MESH_H_

#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
bool ghosted_entity_check(MeshEntity& e, bool throw_error);

//-----------------------------------------------------------------------------
bool shared_entity_check(MeshEntity& e, bool throw_error);

//-----------------------------------------------------------------------------
bool ghosted_entities_check(Mesh& boundary, uint dim, bool throw_error);

//-----------------------------------------------------------------------------
bool shared_entities_check(Mesh& boundary, uint dim, bool throw_error);

//-----------------------------------------------------------------------------
bool interior_boundary_entities_check(Mesh& mesh, uint dim, bool throw_error);

}

#endif /* __PTESTS_MESH_H_ */
