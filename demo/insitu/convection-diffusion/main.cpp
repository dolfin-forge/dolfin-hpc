// Copyright (C) 2010 Jeannette Spuehler.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2017.

#include "ConvectionDiffusion.h"
#include <sstream>
#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

struct Source : public Value<Source,2>
{
  void eval(real *value, const real *x) const
  {
    if(fabs(x[0]-0.75)<0.1 && fabs(x[1]-0.5)<0.1) {
      value[0]= 1.0;
      value[1]= 1.0;
    }
    else{
      value[0]= 0.0;
      value[1]= 0.0;
    }
  }
};

struct FluidVelocity : public Value<FluidVelocity,2>
{
  void eval(real *value, const real *x) const
  {
    value[0]= 1.0*(x[1]-0.5);
    value[1]= 1.0*(-x[0]+0.5);    
  }
};

//-----------------------------------------------------------------------------

struct DirichletFunction : public Value<DirichletFunction, 2>
{
  void eval(real *value, const real *x) const
  {
    value[0] = 0.0;
    value[1] = 0.0;
  }
};

struct DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return on_boundary;
  }
};

//-----------------------------------------------------------------------------

// Pipeline to render the solution
class RenderU : public libsimPipeline
{
public: 

  void exec(real t, uint step) const
  {
    
    std::stringstream filename;
    filename << "solution" << std::setfill('0') << std::setw(6) 
	     << step << ".png" << std::ends;
    
    VisItAddPlot("Pseudocolor", "U");
    VisItDrawPlots();
    VisItSaveWindow(filename.str().c_str(), 800, 600, VISIT_IMAGEFORMAT_PNG);
    VisItDeleteActivePlots();
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
  
  libsimInterface::initBatch();

  // Create mesh
  Mesh mesh("UnitSquareMesh_32x32.xml");
  mesh.refine();
  mesh.refine();

  // Expose to mesh to libsim
  libsimInterface::addData(mesh);

  Analytic<Source> f(mesh);
  Analytic<FluidVelocity> W(mesh);
  Constant alpha(0.0);
  Constant nu(0.0001);
  Constant s(0.0158);
  Constant k(0.02);		

  // Add pipeline to the insitu interface
  RenderU render_u;
  libsimInterface::addPipeline(render_u); 

  Analytic<DirichletFunction> u0(mesh);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  

  ConvectionDiffusion::BilinearForm a(mesh, W, alpha, nu, s, k);

  Function u(a.trial_space());
  Function up(a.trial_space());

  // Expose the function named U to libsim
  libsimInterface::addData(u, "U", "Mesh");

  ConvectionDiffusion::LinearForm L(mesh, f, W, up, alpha, nu, s, k);


  Matrix A;
  Vector b;
  a.assemble(A, true);
  
  KrylovSolver solver(bicgstab, bjacobi); 
  for (uint step = 0; step < 10; step++) 
  {
    L.assemble(b, step == 0);
    bc.apply(A, b, a);
    solver.solve(A, u.vector(), b);
    u.sync();
    up = u;

    // Render all pipelines
    libsimInterface::batchRender(0.02 * (real) step, step);
  }

  dolfin_finalize();
  return 0;
  
}
  
  
  
  
