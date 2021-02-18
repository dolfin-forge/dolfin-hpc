
#ifndef __DOLFIN_FUNCTION_REAL_H_
#define __DOLFIN_FUNCTION_REAL_H_

#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/ValueSpace.h>
#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

template < size_t I = 1, size_t J = 1 >
class Real : public Coefficient
{
public:
  /// Default constructor
  Real();

  /// Constructor
  explicit Real( real value );

  /// Destructor
  ~Real();

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  auto evaluate( real * values, const real *, const ufc::cell & ) const -> void override;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point in cell
  auto evaluate( size_t, real * values, const real *, const ufc::cell & ) const -> void override;

  /// Evaluate function at given point
  auto eval( real * values, const real * ) const -> void override;

  /// Return the rank of the value space
  auto rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto dim( size_t i ) const -> size_t override;

  /// Value size
  auto value_size() const -> size_t override;

  /// Assign constant real number
  auto operator=( real value ) -> Real &;

  /// Multiply by constant real number
  auto operator+=( real value ) -> Real &;

  /// Add a constant real number
  auto operator-=( real value ) -> Real &;

  /// Substract a constant real number
  auto operator*=( real value ) -> Real &;

  /// Divide by constant real number
  auto operator/=( real value ) -> Real &;

  /// Accessors
  auto operator[]( size_t i ) -> real &;

  ///
  auto operator[]( size_t i ) const -> real const &;

  ///
  auto operator()( Time const & ) const -> Real< I, J > const &;

  /// Interpolate function to finite element space on cell
  auto interpolate( real *                      coefficients,
                    UFCCell const &             cell,
                    ufc::finite_element const & finite_element ) const
    -> void override;

  /// Interpolate function to finite element space on facet
  auto interpolate( real *                      coefficients,
                    UFCCell const &             cell,
                    ufc::finite_element const & finite_element,
                    size_t ) const -> void override;

  /// Synchronize
  auto sync() -> void override;

  /// Display basic information
  auto disp() const -> void override;

private:
  auto sync( Time const & ) -> void override;

  real value_[I * J];
};

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline Real< I, J >::Real()
  : value_()
{
  std::fill_n( value_, I * J, 0.0 );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline Real< I, J >::Real( real value )
  : value_()
{
  ( *this ) = value;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline Real< I, J >::~Real()
{
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::evaluate( real * values,
                                    const real *,
                                    const ufc::cell & ) const -> void
{
  std::copy( value_, value_ + I * J, values );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::evaluate( size_t,
                                    real * values,
                                    const real *,
                                    const ufc::cell & ) const -> void
{
  std::copy( value_, value_ + I * J, values );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::eval( real * values, const real * ) const -> void
{
  std::copy( value_, value_ + I * J, values );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::rank() const -> size_t
{
  return ValueSpace< I, J >::rank();
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::dim( size_t i ) const -> size_t
{
  return ValueSpace< I, J >::dim( i );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::value_size() const -> size_t
{
  return ValueSpace< I, J >::value_size();
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator=( real value ) -> Real &
{
  std::fill_n( value_, I * J, value );
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator+=( real value ) -> Real &
{
  for ( size_t i = 0; i < I * J; ++i )
    value_[i] += value;
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator-=( real value ) -> Real &
{
  for ( size_t i = 0; i < I * J; ++i )
    value_[i] -= value;
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator*=( real value ) -> Real &
{
  for ( size_t i = 0; i < I * J; ++i )
    value_[i] *= value;
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator/=( real value ) -> Real &
{
  for ( size_t i = 0; i < I * J; ++i )
    value_[i] /= value;
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator[]( size_t i ) -> real &
{
  return value_[i];
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator[]( size_t i ) const -> real const &
{
  return value_[i];
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::operator()( Time const & ) const
  -> Real< I, J > const &
{
  return *this;
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto
  Real< I, J >::interpolate( real *                      coefficients,
                             UFCCell const &             cell,
                             ufc::finite_element const & finite_element ) const
  -> void
{
  dolfin_assert( coefficients != nullptr );

  finite_element.evaluate_dofs( coefficients, *this, cell.coordinates.data(),
                                cell.orientation, cell );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::interpolate( real *                      coefficients,
                                       UFCCell const &             cell,
                                       ufc::finite_element const & finite_element,
                                       size_t ) const -> void
{
  dolfin_assert( coefficients != nullptr );

  finite_element.evaluate_dofs( coefficients, *this, cell.coordinates.data(),
                                cell.orientation, cell );
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::sync() -> void
{
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::disp() const -> void
{
  section( "Real" );
  message( "ValueSpace<%u,%u>", I, J );
  size_t k = 0;
  for ( size_t i = 0; i < I; ++i )
  {
    cout << "\n";
    for ( size_t j = 0; j < J; ++j, ++k )
    {
      cout << "\t" << value_[k] << "";
    }
  }
  cout << "\n";
  end();
  skip();
}

//-----------------------------------------------------------------------------

template < size_t I, size_t J >
inline auto Real< I, J >::sync( Time const & ) -> void
{
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_REAL_H_ */
