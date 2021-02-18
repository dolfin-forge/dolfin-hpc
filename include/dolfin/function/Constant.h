
#ifndef __DOLFIN_CONSTANT_H_
#define __DOLFIN_CONSTANT_H_

#include <dolfin/fem/Coefficient.h>
#include <dolfin/log/log.h>

namespace dolfin
{

class Time;

class Constant : public Coefficient
{

public:
  /// Constructor
  Constant();

  /// Constructor
  Constant( real value );

  /// Destructor
  ~Constant();

  /// Assignment
  auto operator=( real const & value ) -> Constant &;

  /// Equality
  auto operator==( Constant const & other ) -> bool;

  /// Multiply by constant real number
  auto operator+=( real const & value ) -> Constant &;

  /// Add a constant real number
  auto operator-=( real const & value ) -> Constant &;

  /// Substract a constant real number
  auto operator*=( real const & value ) -> Constant &;

  /// Divide by constant real number
  auto operator/=( real const & value ) -> Constant &;

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  auto evaluate( real * values, const real *, const ufc::cell & ) const -> void override;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point in cell
  auto evaluate( size_t n, real * values, const real *, const ufc::cell & ) const -> void override;

  /// Evaluate function at given point
  auto eval( real * values, const real * ) const -> void override;

  /// Return the rank of the value space
  auto rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto dim( size_t ) const -> size_t override;

  /// Value size
  auto value_size() const -> size_t override;

  ///
  auto operator()( Time const & ) const -> Constant const &;

  /// Interpolate function to finite element space on cell
  auto interpolate( real * coefficients,
                    UFCCell const &,
                    ufc::finite_element const & finite_element ) const
    -> void override;

  /// Interpolate function to finite element space on facet
  auto interpolate( real * coefficients,
                    UFCCell const &,
                    ufc::finite_element const & finite_element,
                    size_t ) const -> void override;

  /// Synchronize
  auto sync() -> void override;

  /// Display basic information
  auto disp() const -> void override;

private:
  auto sync( Time const & ) -> void override;

  real value_;
};

//-----------------------------------------------------------------------------

inline Constant::Constant()
  : value_( 0.0 )
{
}

//-----------------------------------------------------------------------------

inline Constant::Constant( real value )
  : value_( value )
{
}

//-----------------------------------------------------------------------------

inline Constant::~Constant()
{
}

//-----------------------------------------------------------------------------

inline auto Constant::operator=( real const & value ) -> Constant &
{
  value_ = value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator==( Constant const & other ) -> bool
{
  if ( value_ == other.value_ )
  {
    return true;
  }
  return false;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator+=( real const & value ) -> Constant &
{
  value_ += value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator-=( real const & value ) -> Constant &
{
  value_ -= value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator*=( real const & value ) -> Constant &
{
  value_ *= value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator/=( real const & value ) -> Constant &
{
  value_ /= value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Constant::evaluate( real * values,
                                const real *,
                                const ufc::cell & ) const -> void
{
  values[0] = value_;
}

//-----------------------------------------------------------------------------

inline auto Constant::evaluate( size_t n,
                                real * values,
                                const real *,
                                const ufc::cell & ) const -> void
{
  std::fill_n( values, n, value_ );
}

//-----------------------------------------------------------------------------

inline auto Constant::eval( real * values, const real * ) const -> void
{
  values[0] = value_;
}

//-----------------------------------------------------------------------------

inline auto Constant::rank() const -> size_t
{
  return 0;
}

//-----------------------------------------------------------------------------

inline auto Constant::dim( size_t ) const -> size_t
{
  return 1;
}

//-----------------------------------------------------------------------------

inline auto Constant::value_size() const -> size_t
{
  return 1;
}

//-----------------------------------------------------------------------------

inline auto Constant::operator()( Time const & ) const -> Constant const &
{
  // No-op
  return *this;
}

//-----------------------------------------------------------------------------

inline auto
  Constant::interpolate( real * coefficients,
                         UFCCell const &,
                         ufc::finite_element const & finite_element ) const
  -> void
{
  dolfin_assert( coefficients != nullptr );
  for ( size_t i = 0; i < finite_element.space_dimension(); ++i )
  {
    coefficients[i] = value_;
  }
}

//-----------------------------------------------------------------------------

inline auto Constant::interpolate( real * coefficients,
                                   UFCCell const &,
                                   ufc::finite_element const & finite_element,
                                   size_t ) const -> void
{
  dolfin_assert( coefficients );
  for ( size_t i = 0; i < finite_element.space_dimension(); ++i )
  {
    coefficients[i] = value_;
  }
}

//-----------------------------------------------------------------------------

inline auto Constant::sync() -> void
{
}

//-----------------------------------------------------------------------------

inline auto Constant::disp() const -> void
{
  section( "Constant" );
  message( "Value : %f", value_ );
  end();
  skip();
}

//-----------------------------------------------------------------------------

inline auto Constant::sync( Time const & ) -> void
{
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_CONSTANT_H_ */
