#ifndef __LICORNE_FUNCTION_REAL_H_
#define __LICORNE_FUNCTION_REAL_H_

#include <dolfin/function/Value.h>

namespace dolfin
{

class Real : public Value<Real>
{

public:

  /// Default constructor
  Real() :
    value_(0.0)
  {
  }

  /// Constructor
  explicit Real(real value) :
    value_(value)
  {
  }

  /// Destructor
  ~Real()
  {
  }

  /// Evaluate value
  void eval(real* values, real const* x) const
  {
    values[0] = value_;
  }

  /// Assign constant real number
  inline Real& operator=(real const& value)
  {
    value_ = value;
    return *this;
  }

  /// Multiply by constant real number
  inline Real& operator+=(real const& value)
  {
    value_ += value;
    return *this;
  }

  /// Add a constant real number
  inline Real& operator-=(real const& value)
  {
    value_ -= value;
    return *this;
  }

  /// Substract a constant real number
  inline Real& operator*=(real const& value)
  {
    value_ *= value;
    return *this;
  }

  /// Divide by constant real number
  inline Real& operator/=(real const& value)
  {
    value_ /= value;
    return *this;
  }

  /// Cast to real number
  inline operator real() const
  {
    return value_;
  }

  ///
  inline void str()
  {
    message("Real: %f", value_);
  }

private:

  real value_;

};

} /* namespace licorne */

#endif /* __LICORNE_FUNCTION_REAL_H_ */
