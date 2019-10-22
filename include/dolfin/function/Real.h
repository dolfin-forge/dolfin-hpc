#ifndef __LICORNE_FUNCTION_REAL_H_
#define __LICORNE_FUNCTION_REAL_H_

#include <dolfin/fem/Coefficient.h>
#include <dolfin/function/ValueSpace.h>

namespace dolfin
{

template<uint I = 1, uint J = 1>
class Real : public Coefficient
{

  static const ValueSpace<I, J> VS_;

public:

  /// Default constructor
  Real() :
    value_()
  {
    std::fill_n(value_, I * J, 0.0);
  }

  /// Constructor
  explicit Real(real value) :
    value_()
  {
    (*this) = value;
  }

  /// Destructor
  ~Real()
  {
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate(real* values, const real*, const ufc::cell&) const
  {
    std::copy(value_, value_ + I * J, values);
  }

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate(uint, real* values, const real*, const ufc::cell&) const
  {
    std::copy(value_, value_ + I * J, values);
  }

  /// Evaluate function at given point
  inline void eval(real* values, const real*) const
  {
    std::copy(value_, value_ + I * J, values);
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

  /// Value size
  inline uint value_size() const
  {
    return VS_.value_size();
  }

  /// Assign constant real number
  inline Real& operator=(real value)
  {
    std::fill_n(value_, I * J, value);
    return *this;
  }

  /// Multiply by constant real number
  inline Real& operator+=(real value)
  {
    for(uint i = 0; i < I * J; ++i)
      value_[i] += value;
    return *this;
  }

  /// Add a constant real number
  inline Real& operator-=(real value)
  {
    for(uint i = 0; i < I * J; ++i)
      value_[i] -= value;
    return *this;
  }

  /// Substract a constant real number
  inline Real& operator*=(real value)
  {
    for(uint i = 0; i < I * J; ++i)
      value_[i] *= value;
    return *this;
  }

  /// Divide by constant real number
  inline Real& operator/=(real value)
  {
    for(uint i = 0; i < I * J; ++i)
      value_[i] /= value;
    return *this;
  }

  /// Accessors
  real&       operator[](uint i)       { return value_[i]; }
  real const& operator[](uint i) const { return value_[i]; }

  ///
  inline Real<I, J> const& operator()(Time const&) const
  {
    // No-op
    return *this;
  }

  /// Interpolate function to finite element space on cell
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell&) const
  {
    finite_element.evaluate_dofs(coefficients, *this, cell);
  }

  /// Interpolate function to finite element space on facet
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell&, uint) const
  {
    finite_element.evaluate_dofs(coefficients, *this, cell);
  }

  /// Synchronize
  inline void sync()
  {
    // Do nothing
  }

  /// Display basic information
  inline void disp() const
  {
    section("Real");
    message("ValueSpace<%u,%u>", I, J);
    uint k = 0;
    for (uint i = 0; i < I; ++i)
    {
      cout << "\n";
      for (uint j = 0; j < J; ++j, ++k)
      {
        cout << "\t" << value_[k] << "";
      }
    }
    cout << "\n";
    end();
    skip();
  }

private:

  void sync(Time const&) { /* No-op */ }

  real value_[I * J];

};

} /* namespace licorne */

#endif /* __LICORNE_FUNCTION_REAL_H_ */
