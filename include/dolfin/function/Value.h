
#ifndef __DOLFIN_FUNCTION_VALUE_H_
#define __DOLFIN_FUNCTION_VALUE_H_

#include <dolfin/evolution/TimeDependent.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/ValueSpace.h>
#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

template < class T, size_t I = 1, size_t J = 1 >
class Value : public Expression, public TimeDependent
{

public:
  ///
  Value()
    : Expression()
    , TimeDependent()
  {
  }

  /// Evaluate expression at given point
  inline auto eval( real * values, real const * x ) const -> void override
  {
    static_cast< T const * >( this )->eval( values, x );
  }

  /// Return the rank of the value space
  inline auto rank() const -> size_t override
  {
    return ValueSpace< I, J >::rank();
  }

  /// Return the dimension of the value space for axis i
  inline auto dim( size_t i ) const -> size_t override
  {
    return ValueSpace< I, J >::dim( i );
  }

  /// Return value size (allow overloading to avoid recomputation)
  inline auto value_size() const -> size_t override
  {
    return ValueSpace< I, J >::value_size();
  }

  //---------------------------------------------------------------------------

  /// Value implements the time dependency
  inline Value< T, I, J > const & operator()( Time const & t ) const
  {
    TimeDependent::operator()( t );
    return *this;
  }

  //---------------------------------------------------------------------------

  ///
  virtual void disp() const
  {
    section( "Value" );
    message( "rank       : %d", this->rank() );
    std::stringstream ss;
    ss << "(";
    for ( size_t i = 0; i <= this->rank(); ++i )
    {
      ss << this->dim( i ) << ", ";
    }
    ss << ")";
    message( "dimensions : %s", ss.str().c_str() );
    message( "value size : %d", this->value_size() );
    end();
    skip();
  }
};

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_VALUE_H_ */
