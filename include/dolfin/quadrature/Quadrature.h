// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_QUADRATURE_H
#define __DOLFIN_QUADRATURE_H

#include <dolfin/common/types.h>

namespace dolfin
{
  
  class Quadrature
  {
  public:
    
    /// Constructor
    Quadrature(unsigned int n);

    /// Destructor
    virtual ~Quadrature();
    
    /// Return number of quadrature points
    auto size() const -> unsigned int;

    /// Return quadrature point
    auto point(unsigned int i) const -> real;

    /// Return quadrature weight
    auto weight(unsigned int i) const -> real;

    /// Return sum of weights (length, area, volume)
    auto measure() const -> real;
    
    /// Display quadrature data
    virtual void disp() const = 0;

  protected:
    
    uint n;        // Number of quadrature points
    real* points;  // Quadrature points
    real* weights; // Quadrature weights
    real m;        // Sum of weights
    
  };
  
}

#endif
