// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2007-08-20
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 500*exp(-((x-0.5)^2 + (y-0.5)^2)/0.02)
//
// and boundary conditions given by
//
//     u(x, y)     = 0               for x = 0
//     du/dn(x, y) = 25 sin(5 pi y)  for x = 1
//     du/dn(x, y) = 0               otherwise
//
// The results is rendered insitu using libsim

#include "Poisson.h"

#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

// Source term
struct Source : public Value<Source>
{
  void eval(real * value, const real* x) const
  {
    real dx = x[0] - 0.5;
    real dy = x[1] - 0.5;
    value[0] = 500.0*std::exp(-(dx*dx + dy*dy)/0.02);
  }
};

// Neumann boundary condition
struct Flux : public Value<Flux>
{
  void eval(real * value, const real* x) const
  {
    if (x[0] > DOLFIN_EPS)
      value[0] = 25.0*std::sin(5.0*DOLFIN_PI*x[1]);
    else
      value[0] = 0.0;
  }
};

//-----------------------------------------------------------------------------

// Pipeline to render the mesh
class RenderMesh : public libsimPipeline
{
public:

  void exec(real t, uint step) const
  {
    VisItAddPlot("Mesh", "Mesh");
    VisItDrawPlots();
    VisItSaveWindow("./mesh.png", 800, 600, VISIT_IMAGEFORMAT_PNG);
    VisItDeleteActivePlots();
  }
};

// Pipeline to render the solution
class RenderU : public libsimPipeline
{
public:

  void exec(real t, uint step) const
  {
    VisItAddPlot("Pseudocolor", "U");
    VisItDrawPlots();
    VisItSaveWindow("./solution.png", 800, 600, VISIT_IMAGEFORMAT_PNG);
    VisItDeleteActivePlots();
  }
};

// Pipline to render both mesh and solution in the same plot
class RenderAll : public libsimPipeline
{
public:

  void exec(real t, uint step) const
  {
    VisItAddPlot("Mesh", "Mesh");
    VisItAddPlot("Pseudocolor", "U");
    VisItDrawPlots();
    VisItSaveWindow("./mesh_and_solution.png",
		    800, 600, VISIT_IMAGEFORMAT_PNG);
    VisItDeleteActivePlots();
  }
};

// Pipline to render a predefined session
class RenderSession : public libsimPipeline
{
public:

  void exec(real t, uint step) const
  {
    VisItRestoreSession("./poisson.session");
    VisItSaveWindow("./session.png",
		    800, 600, VISIT_IMAGEFORMAT_PNG);
  }
};

//-----------------------------------------------------------------------------

// Sub domain for Dirichlet boundary condition
struct DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS && on_boundary;
  }
};

//-----------------------------------------------------------------------------

int main(int argc, char **argv)
{

  dolfin_init(argc, argv);

  // VisIt directory should point to the top installation dir.
  // dolfin_set("VisIt directory","/opt/visit/2.10.0");

  const std::string visit_path = dolfin_get("VisIt directory");
  if (visit_path == "")
    error("Parameter 'VisIt directory' not set properly");

  libsimInterface::init(libsimInterface::batch);

  // Create mesh
  Mesh mesh("UnitSquareMesh_32x32.bin");

  // Expose a mesh named Mesh to libsim
  libsimInterface::addData(mesh, "Mesh");

  // Create coefficients
  Analytic<Source>  f(mesh);
  Analytic<Flux>    g(mesh);
  RenderU render_u;
  RenderMesh render_mesh;
  RenderAll render_all;
  RenderSession render_session;

  // Add all pipelines to the insitu interface
  libsimInterface::addPipeline(render_u);
  libsimInterface::addPipeline(render_mesh);
  libsimInterface::addPipeline(render_all);
  libsimInterface::addPipeline(render_session);

  // Create boundary condition
  Constant u0(0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);

  // Define PDE
  Poisson::BilinearForm a(mesh);
  Poisson::LinearForm   L(mesh, f, g);

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble(A, true);
  L.assemble(b, true);
  bc.apply(A, b, a);


  Function u(a.trial_space());

  // Expose the function named U to libsim
  libsimInterface::addData(u, "U", "Mesh");

  KrylovSolver solver(bicgstab, bjacobi);
  solver.solve(A, u.vector(), b);
  u.sync();

  // Render all pipelines
  libsimInterface::batchRender();

  return 0;
}
