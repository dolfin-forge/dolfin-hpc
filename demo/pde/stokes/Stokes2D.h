// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-09
// Last changed: 2007-07-11
//
// This demo solves the Stokes equations, using stabilized
// first order elements for the velocity and pressure. The
// sub domains for the different boundary conditions used
// in this simulation are computed by the demo program in
// src/demo/mesh/subdomains.

#include <dolfin.h>

using namespace dolfin;

class LeftBoundary : public SubDomain
{
public:

	bool inside(const real* p, bool on_boundary) const
	{
		return on_boundary&& (p[0] <= DOLFIN_EPS);
	}
};

class RightBoundary : public SubDomain
{
public:

	bool inside(const real* p, bool on_boundary) const
	{
		return on_boundary&& (p[0] >= (1 -  DOLFIN_EPS));
	}
};

class NoSlipBoundary : public SubDomain
{
public:

	bool inside(const real* p, bool on_boundary) const
	{	
		return on_boundary&&(p[0] > DOLFIN_EPS)&&(p[0] < (1-DOLFIN_EPS));
	}
};

// Function for no-slip boundary condition for velocity
class NoSlip : public VectorExpression
{
public:

	NoSlip() : VectorExpression(2) {}

	void eval(real* values, const real* x) const
	{
		values[0] = 0.0;
		values[1] = 0.0;
	}
};

// Function for inflow boundary condition for velocity
class Inflow : public VectorExpression
{
public:

	Inflow() : VectorExpression(2) {}

	void eval(real* values, const real* x) const
	{
		values[0] = -sin(x[1]*DOLFIN_PI);
		values[1] = 0.0;
	}
};

