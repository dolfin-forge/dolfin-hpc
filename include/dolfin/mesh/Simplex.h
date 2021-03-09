// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#ifndef __DOLFIN_MESH_SIMPLEX_H
#define __DOLFIN_MESH_SIMPLEX_H

#include <dolfin/mesh/Mesh.h>

#include <dolfin/mesh/celltypes/IntervalCell.h>
#include <dolfin/mesh/celltypes/TetrahedronCell.h>
#include <dolfin/mesh/celltypes/TriangleCell.h>
#include <dolfin/mesh/EuclideanSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
template < size_t D >
struct Simplex : public Mesh
{

private:
	Simplex()
	{
	}
};

//-----------------------------------------------------------------------------
template <>
struct Simplex< 1 > : public Mesh
{
	Simplex< 1 >()
		: Mesh( IntervalCell(), EuclideanSpace( 1 ) )
	{
		IntervalCell C;
		C.create_reference_cell( *this );
	}
};

//-----------------------------------------------------------------------------
template <>
struct Simplex< 2 > : public Mesh
{
	Simplex< 2 >()
		: Mesh( TriangleCell(), EuclideanSpace( 2 ) )
	{
		TriangleCell C;
		C.create_reference_cell( *this );
	}
};

//-----------------------------------------------------------------------------
template <>
struct Simplex< 3 > : public Mesh
{
	Simplex< 3 >()
		: Mesh( TetrahedronCell(), EuclideanSpace( 3 ) )
	{
		TetrahedronCell C;
		C.create_reference_cell( *this );
	}
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SIMPLEX_H */
