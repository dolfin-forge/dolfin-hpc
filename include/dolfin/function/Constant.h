#ifndef __DOLFIN_CONSTANT_H_
#define __DOLFIN_CONSTANT_H_

#include <dolfin/fem/Coefficient.h>

#include <dolfin/log/log.h>

namespace dolfin
{

class Constant : public Coefficient
{

public:

  /// Constructor
  Constant() :
    value_(0.0)
  {
  }

  /// Constructor
  Constant(real value) :
    value_(value)
  {
  }

  /// Destructor
  ~Constant()
  {
  }

  /// Assignment
  Constant& operator=(real const& value)
  {
    value_ = value;
    return *this;
  }

  /// Equality
  bool operator==(Constant const& other)
  {
    if(value_ == other.value_)
    {
      return true;
    }
    return false;
  }

  /// Multiply by constant real number
  inline Constant& operator+=(real const& value)
  {
    value_ += value;
    return *this;
  }

  /// Add a constant real number
  inline Constant& operator-=(real const& value)
  {
    value_ -= value;
    return *this;
  }

  /// Substract a constant real number
  inline Constant& operator*=(real const& value)
  {
    value_ *= value;
    return *this;
  }

  /// Divide by constant real number
  inline Constant& operator/=(real const& value)
  {
    value_ /= value;
    return *this;
  }

  ///
  operator real() const
  {
    return value_;
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate(real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    values[0] = value_;
  }

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate(uint n, real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    std::fill_n(values, n, value_);
  }

  /// Evaluate function at given point
  inline void eval(real* values, const real* x) const
  {
    values[0] = value_;
  }

  /// Return the rank of the value space
  inline uint rank() const
  {
    return 0;
  }

  /// Return the dimension of the value space for axis i
  inline uint dim(uint i) const
  {
    return 1;
  }

  /// Value size
  inline uint value_size() const
  {
    return 1;
  }

  /// Interpolate function to finite element space on cell
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell& dolfin_cell) const
  {
    dolfin_assert(coefficients);
    for (uint i = 0; i < finite_element.space_dimension(); ++i)
    {
      coefficients[i] = value_;
    }
  }

  /// Interpolate function to finite element space on facet
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell& dolfin_cell, uint facet) const
  {
    dolfin_assert(coefficients);
    for (uint i = 0; i < finite_element.space_dimension(); ++i)
    {
      coefficients[i] = value_;
    }
  }

  /// Synchronize
  inline void sync()
  {
    // Do nothing
  }

  /// Display basic information
  inline void disp() const
  {
    section("Constant");
    message("Value : %f", value_);
    end();
    skip();
  }

private:

  real value_;

};

} /* namespace licorne */

#endif /* __DOLFIN_CONSTANT_H_ */
