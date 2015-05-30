#include <dolfin.h>

namespace dolfin
{

// Source term
class Source : public Function
{
public:

  Source(Mesh& mesh) :
      Function(mesh)
  {
  }

  void eval(real * value, const real* x) const
  {
    real dx = x[0] - 0.5;
    real dy = x[1] - 0.5;
    value[0] = 500.0 * exp(-(dx * dx + dy * dy) / 0.02);
  }

  uint rank() const
  {
    return 0;
  }

  uint dim(uint i) const
  {
    return 1;
  }

};

// Neumann boundary condition
class Flux : public Function
{
public:

  Flux(Mesh& mesh) :
      Function(mesh)
  {
  }

  void eval(real * value, const real* x) const
  {
    if (x[0] > DOLFIN_EPS) value[0] = 25.0 * sin(5.0 * DOLFIN_PI * x[1]);
    else value[0] = 0.0;
  }

  uint rank() const
  {
    return 0;
  }

  uint dim(uint i) const
  {
    return 1;
  }

};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS && on_boundary;
  }
};

}
