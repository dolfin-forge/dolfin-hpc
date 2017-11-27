// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.
//
// First added:
// Last changed:

#ifndef __DOLFIN_MESH_SIMPLEX_H
#define __DOLFIN_MESH_SIMPLEX_H

#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>

namespace dolfin
{

template<uint D> struct Simplex;

//-----------------------------------------------------------------------------
template<>
struct Simplex<1> : public Mesh
{
  Simplex<1>() { IntervalCell C; C.create_reference_cell(*this); }
};

//-----------------------------------------------------------------------------
template<>
struct Simplex<2> : public Mesh
{
  Simplex<2>() { TriangleCell C; C.create_reference_cell(*this); }
};

//-----------------------------------------------------------------------------
template<>
struct Simplex<3> : public Mesh
{
  Simplex<3>() { TetrahedronCell C; C.create_reference_cell(*this); }
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SIMPLEX_H */
