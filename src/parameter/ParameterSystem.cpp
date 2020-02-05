// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/parameter/ParameterSystem.h>
#include <dolfin/common/constants.h>

using namespace dolfin;

// Initialize the global parameter database
ParameterSystem ParameterSystem::parameters;

//-----------------------------------------------------------------------------
ParameterSystem::ParameterSystem()
  : ParameterList()
{
#include <dolfin/config/dolfin_config.h>


//--- Linear algebra ---
#ifdef HAVE_PETSC
  add( "linear algebra backend", "PETSc" );
#elif HAVE_JANPACK
  add( "linear algebra backend", "JANPACK" );
#else
  add( "linear algebra backend", "none" );
#endif

  //--- JIT compiler ---
  add( "optimize form",
       false ); // Use optimization -O2 when compiling generated code
  add( "optimize use dof map cache",
       false ); // Store dof maps in cache for reuse
  add( "optimize use tensor cache", false ); // Store tensors in cache for reuse
  add( "optimize", false );                  // All of the above

  //--- General parameters ---

  add( "solution file name", "solution.pvd" );

  //--- Parameters for input/output ---

  add( "output_format", "vtk" );

  //--- Parameters for Krylov solvers ---

  add( "Krylov relative tolerance", 1e-10 );
  add( "Krylov absolute tolerance", 1e-20 );
  add( "Krylov divergence limit", 1e4 );
  add( "Krylov maximum iterations", 10000 );
  add( "Krylov GMRES restart", 30 );
  add( "Krylov shift nonzero", 0.0 );
  add( "Krylov report", true );
  add( "Krylov keep PC", false );
  add( "Krylov error on nonconvergence", true );

  //--- Parameters for AMG solvers ---

  add( "AMG relative tolerance", 1e-10 );
  add( "AMG absolute tolerance", 1e-20 );
  add( "AMG maximum iterations", 10000 );
  add( "AMG pre-smoothing steps", 5 );
  add( "AMG post-smoothing steps", 5 );
  add( "AMG theta", 0.15 );
  add( "AMG levels", 20 );
  add( "AMG keep levels", false );

  //--- Parameter for direct (LU) solver ---
  add( "LU report", true );

  //--- Parameter for PDE solver ---
  add( "PDE linear solver", "iterative" );

  //--- Mesh partitioning ---
  add( "report edge cut", false ); // ?
  add( "Mesh read in serial", false );
#if HAVE_PARMETIS
  add( "Mesh partitioner", "parmetis" );
#elif HAVE_ZOLTAN
  add( "Mesh partitioner", "zoltan" );
#endif

  //--- Load balancing ---
  add( "Load balancer report", false );
  add( "Load balancer redistribute", true );

  //--- Parameters for intersection detection ---
  add( "GTS Tolerance", 0.0 ); // Tolerance of GTS BB
  // define size for trianlge tolerance ("is the point within this triangle?")
  add( "Geometrical Tolerance Interval", 0.0 );
  add( "Geometrical Tolerance Triangle", 0.0 );
  add( "Geometrical Tolerance Tetrahedron", 0.0 );
  add( "SubDomain Geometrical Tolerance", 1e-6 );
  add( "SubDomain Intersect Boundary", true );

  //--- Mesh smoothing ---
  add( "Mesh smoothing restricted by rmin", true );

  //--- Insitu ---
  add( "VisIt directory", "" );

  //--- Insitu ---
  add( "NodeNormal alpha", DOLFIN_PI / 2. );
}
//-----------------------------------------------------------------------------
