// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_COEFFICIENT_H
#define __DOLFIN_COEFFICIENT_H

#include <dolfin/common/types.h>
#include <dolfin/evolution/Time.h>

#include <ufc.h>

namespace dolfin
{

class Mesh;
class Cell;

/// This class serves as a base class/interface for implementations
/// of specific function representations.

class Coefficient : public ufc::function
{
public:

  /// Constructor
  Coefficient()
  {
  }

  /// Destructor
  virtual ~Coefficient()
  {
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  virtual void evaluate(real* values, const real* coordinates,
                        const ufc::cell& cell) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given points in cell
  virtual void evaluate(uint n, real* values, const real* coordinates,
                        const ufc::cell& cell) const = 0;

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

  /// Return the rank of the value space
  virtual uint rank() const = 0;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const = 0;

  // Return the value size
  virtual uint value_size() const = 0;

  /// Interpolate function to finite element space on cell
  virtual void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell) const = 0;

  /// Interpolate function to finite element space on facet
  virtual void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const = 0;

  /// Synchronize values
  virtual void sync() = 0;

  /// Display basic information
  virtual void disp() const = 0;

  //---------------------------------------------------------------------------

  /// Delegate time dependence
  Coefficient& operator()(Time const& t)
  {
    this->sync(t);
    return *this;
  }

  ///
  template<class T>
  bool compatible(T const& other)
  {
    return compatible(*this, other);
  }

  ///
  template<class T>
  static bool compatible(Coefficient const& a, T const& b)
  {
    // Check value shape compatibility
    if (a.rank() != b.rank())
    {
      warning("Coefficient: rank mismatch %u != %u", a.rank(), b.rank());
      return false;
    }
    for (uint i = 0; i < a.rank(); ++i)
    {
      if (a.dim(i) != b.dim(i))
      {
        warning("Coefficient: dim(%u) mismatch %u != %u", i, a.dim(i), b.dim(i));
        return false;
      }
    }
    return true;
  }

private:

  /// Time dependency hook
  virtual void sync(Time const& t) = 0;

};

} /* namespace dolfin */

#endif /* __DOLFIN_COEFFICIENT_H */
