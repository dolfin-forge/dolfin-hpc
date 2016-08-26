// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_UFC_COEFFICENT_H
#define __DOLFIN_UFC_COEFFICENT_H

#include <dolfin/fem/Coefficient.h>

#include <dolfin/fem/UFCCell.h>
#include <dolfin/log/log.h>

namespace dolfin
{

/**
 *  @class  UFCCoefficient
 *
 *  @brief  This class serves as a base class/interface for implementations of
 *          user-defined coefficients using UFC cell and facet information.
 */

template<class T, class ValueSpaceImpl>
class UFCCoefficient : public Coefficient
{

  static const ValueSpaceImpl VS_;

public:

  //--- UFC INTERFACE ---------------------------------------------------------

  ///
  inline void evaluate(real* values, const real* coordinates, const ufc::cell& cell) const
  {
    static_cast<T const *>(this)->evaluate(values, coordinates, static_cast<UFCCell const&>(cell));
  }

  //--- Coefficient INTERFACE -------------------------------------------------

  ///
  inline void evaluate(uint n, real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    for (uint i = 0; i < n; ++i)
    {
      static_cast<T const *>(this)->evaluate(&values[i*value_size()],
                                             &coordinates[i*cell.geometric_dimension],
                                             static_cast<UFCCell const&>(cell));
    }
  }

  /// Evaluate function at given point coordinate, cell is searched for
  void eval(real* values, const real* x) const
  {
    error("UFCCoefficient : eval unimplemented");
  }

  /// Return the rank of the value space
  inline uint rank() const
  {
    return VS_.rank();
  }

  /// Return the dimension of the value space for axis i
  inline uint dim(uint i) const
  {
    return VS_.dim(i);
  }

  // Return the value size
  inline uint value_size() const
  {
    return VS_.value_size();
  }

  /// Interpolate function to finite element space on cell
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell& dolfin_cell) const
  {
    dolfin_assert(coefficients);
    finite_element.evaluate_dofs(coefficients, *this, cell);
  }

  /// Interpolate function to finite element space on facet
  inline  void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const
  {
    this->facet_ = facet;
    dolfin_assert(coefficients);
    finite_element.evaluate_dofs(coefficients, *this, cell);
    this->facet_ = -1;
  }

  /// Display basic information
  void disp() const
  {
    section("UFCCoefficient");
    section("ValueSpace");
    message("rank       : %u", this->rank());
    message("value size : %u", this->value_size());
    endblock();
    endblock();
  }

  /// Synchronize
  void sync() { /* No-op */ }

protected:

  inline int facet() const { return facet_; }

private:

  /// Default time dependency hook
  virtual void sync(Time const& t) { /* No-op */ }

  mutable int facet_;

};

//-----------------------------------------------------------------------------
template<class T, class ValueSpaceImpl>
ValueSpaceImpl const UFCCoefficient<T, ValueSpaceImpl>::VS_ = ValueSpaceImpl();

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_COEFFICENT_H */
