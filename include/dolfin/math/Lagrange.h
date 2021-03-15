// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LAGRANGE_H
#define __DOLFIN_LAGRANGE_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>

namespace dolfin
{
  /// Lagrange polynomial (basis) with given degree q determined by n = q + 1 nodal points.
  ///
  /// Example: q = 1 (n = 2)
  ///
  ///   Lagrange p(1);
  ///   p.set(0, 0.0);
  ///   p.set(1, 1.0);
  ///
  /// This creates a Lagrange polynomial (actually two Lagrange polynomials):
  ///
  ///   p(0,x) = 1 - x   (one at x = 0, zero at x = 1)
  ///   p(1,x) = x       (zero at x = 0, one at x = 1)

  class Lagrange
  {
  public:
 
    /// Constructor
    Lagrange(unsigned int q);
    
    /// Copy constructor
    Lagrange(const Lagrange& p);
    
    /// Destructor
    ~Lagrange();
    
    /// Specify point
    void set(unsigned int i, real x);

    /// Return number of points
    auto size() const -> unsigned int;
    
    /// Return degree
    auto degree() const -> unsigned int;
    
    /// Return point
    auto point(unsigned int i) const -> real;
    
    /// Return value of polynomial i at given point x
    auto operator() (unsigned int i, real x) -> real;
    
    /// Return value of polynomial i at given point x
    auto eval(unsigned int i, real x) -> real;

    /// Return derivate of polynomial i at given point x
    auto ddx(unsigned int i, real x) -> real;

    /// Return derivative q (a constant) of polynomial
    auto dqdx(unsigned int i) -> real;

    /// Output
    friend auto operator<<(LogStream& stream, const Lagrange& p) -> LogStream&;
    void disp() const;
    
  private:

    void init();
    
    unsigned int q;
    unsigned int n;
    real* points; 
    real* constants;

  };

  auto operator<<(LogStream& stream, const Lagrange& p) -> LogStream&;

}

#endif
