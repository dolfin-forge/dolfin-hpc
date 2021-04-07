#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/la/Matrix.h>

#define MAT_SIZE 10

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_init_mat )
{
  Matrix A;
  A.init(MAT_SIZE, MAT_SIZE);

  ck_assert(A.size(0) == MAT_SIZE);
  ck_assert(A.size(1) == MAT_SIZE);

}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
#ifdef HAVE_TRILINOS
#include <dolfin/la/trilinos/TrilinosMatrix.h>
#undef fail
#include <dolfin/la/trilinos/TrilinosKrylovSolver.h>
#include "../../demo/pde/poisson/Poisson.h"
#include <dolfin/function/Function.h>
#include <dolfin/function/Analytic.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Value.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/mesh/unitmeshes/UnitSquare.h>
#include <dolfin/common/constants.h>

// example from demo/pde/poisson

// Source term
struct Source : public Value< Source >
{
  void eval( real * value, const real * x ) const
  {
    real dx  = x[0] - 0.5;
    real dy  = x[1] - 0.5;
    value[0] = 500.0 * std::exp( -( dx * dx + dy * dy ) / 0.02 );
  }
};

// Neumann boundary condition
struct Flux : public Value< Flux >
{
  void eval( real * value, const real * x ) const
  {
    if ( x[0] > DOLFIN_EPS )
      value[0] = 25.0 * std::sin( 5.0 * DOLFIN_PI * x[1] );
    else
      value[0] = 0.0;
  }
};

//------------------------------------------------------------------------

// Sub domain for Dirichlet boundary condition
struct DirichletBoundary : public SubDomain
{
  bool inside( const real * x, bool on_boundary ) const
  {
    return x[0] < DOLFIN_EPS && on_boundary;
  }
};

DOLFIN_START_TEST( test_trilinos_mat )
{
  trilinos::Matrix x1;
  trilinos::Matrix x2( MAT_SIZE, MAT_SIZE );
  trilinos::Matrix x3( x2 );

  x1 = x3;

  ck_assert(x1.size(0) == MAT_SIZE);
  ck_assert(x1.size(1) == MAT_SIZE);

  x1.zero();
  x1.apply();
  x1.disp();

  {
    //------------------------------------------------------------------------

  	dolfin_set( "linear algebra backend", "Trilinos" );

    // Create mesh
    // Mesh mesh("./demo/pde/poisson/UnitSquareMesh_32x32.bin");
    UnitSquare mesh( 32, 32 );

    // Create coefficients
    Analytic< Source > f( mesh );
    Analytic< Flux >   g( mesh );

    // Create boundary condition
    Constant          u0( 0.0 );
    DirichletBoundary boundary;
    DirichletBC       bc( u0, mesh, boundary );

    // Define PDE
    Poisson::BilinearForm a( mesh );
    Poisson::LinearForm   L( mesh, f, g );

    // Solve PDE
    Matrix A;
    Vector b;

    a.assemble( A, true );
    L.assemble( b, true );
    bc.apply( A, b, a );

    Function u( a.trial_space() );
    u.vector().disp();
    // KrylovSolver solver( bicgstab, bjacobi );

    // solver.solve( A, u.vector(), b );
    // u.sync();

    // BinaryFile( "u.bin" ) << u;

    message( "vector l2  norm: %e", u.vector().norm() );
    message( "vector inf norm: %e", u.vector().max() );
  }
  {
    message( "solve test" );
    dolfin_set( "linear algebra backend", "Trilinos" );

    trilinos::Matrix A( MAT_SIZE, MAT_SIZE );
    trilinos::Vector x( MAT_SIZE );
    trilinos::Vector b( MAT_SIZE );

    A.zero();
    A.apply();
    x.zero();
    x.apply();
    b.zero();
    b.apply();

    trilinos::KrylovSolver solver1( SolverType::cg, PreconditionerType::relax );
    solver1.disp();
    solver1.solve( A, x, b );

    trilinos::KrylovSolver solver2( SolverType::gmres, PreconditionerType::cheb );
    solver2.disp();
    solver2.solve( A, x, b );

    trilinos::KrylovSolver solver3( SolverType::bicgstab, PreconditionerType::ilut );
    solver3.disp();
    solver3.solve( A, x, b );

  }
}DOLFIN_END_TEST
#else
DOLFIN_START_TEST( test_trilinos_mat )
{
}DOLFIN_END_TEST
#endif
//-----------------------------------------------------------------------------

#endif
