#ifndef __LICORNE_FUNCTION_VALUE_H_
#define __LICORNE_FUNCTION_VALUE_H_

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include <dolfin/evolution/Time.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/ValueSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

template<class T, uint I = 1, uint J = 1>
class Value : public Expression
{
  static ValueSpace<I, J> const VS_;

public:

  ///
  Value() :
      Expression(),
      t_(0.0)
  {
  }

  /// Evaluate expression at given point
  void eval(real* values, real const* x) const
  {
    static_cast<T const *>(this)->eval(values, x);
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

  /// Return value size (allow overloading to avoid recomputation)
  inline uint value_size() const
  {
    return VS_.value_size();
  }

  ///
  inline Value<T, I,J> const& operator()(Time const& t) const
  {
    t_ = t.clock();
    return *this;
  }

  ///
  inline Value<T, I,J> const& operator()(Time const& t, real* values, real const* x) const
  {
    t_ = t.clock();
    eval(values, x);
    return *this;
  }

  ///
  inline real time() const { return t_; }

  ///
  virtual void disp() const
  {
    section("Value");
    message("rank       : %d", this->rank());
    std::stringstream ss;
    ss << "(";
    for (uint i = 0; i <= this->rank(); ++i)
    {
      ss << this->dim(i) << ", ";
    }
    ss << ")";
    message("dimensions : %s", ss.str().c_str());
    message("value size : %d", this->value_size());
    end();
    skip();
  }

private:

  mutable real t_;

};

//
template<class T, uint I, uint J>
ValueSpace<I, J> const Value<T, I,J>::VS_ = ValueSpace<I, J>();

} /* namespace licorne */

#endif /* __LICORNE_FUNCTION_VALUE_H_ */
