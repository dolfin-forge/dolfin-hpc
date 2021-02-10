// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_COEFFICENT_H
#define __DOLFIN_UFC_COEFFICENT_H

#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/ValueSpace.h>
#include <dolfin/log/log.h>

namespace dolfin
{

/**
 *  @class  UserCoefficient
 *
 *  @brief  This class serves as a base class/interface for implementations of
 *          user-defined coefficients using UFC cell and facet information.
 */

template < typename Functor, size_t I = 1, size_t J = 1 >
class UserCoefficient : public Coefficient
{

public:

  //--- UFC INTERFACE ---------------------------------------------------------

  ///
  auto evaluate( real *            values,
                 real const *      coordinates,
                 ufc::cell const & cell ) const -> void override;

  //--- Coefficient INTERFACE -------------------------------------------------

  ///
  auto evaluate( size_t            n,
                 real *            values,
                 real const *      coordinates,
                 ufc::cell const & cell ) const -> void override;

  /// Evaluate function at given point coordinate, cell is searched for
  auto eval( real *, const real * ) const -> void override;

  /// Return the rank of the value space
  auto rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto dim( size_t i ) const -> size_t override;

  /// Return value size (allow overloading to avoid recomputation)
  auto value_size() const -> size_t override;

  /// Interpolate function to finite element space on cell
  auto interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell & ) const -> void override;

  /// Interpolate function to finite element space on facet
  auto interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell &,
                    size_t ) const -> void override;

  /// Display basic information
  auto disp() const -> void override;

  /// Synchronize
  auto sync() -> void override;

protected:
  /// Default time dependency hook
  virtual auto sync( Time const & ) -> void override;

  mutable UFCCell ufc_cell_;
};

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto
  UserCoefficient< Functor, I, J >::evaluate( real *            values,
                                              const real *      coordinates,
                                              const ufc::cell & cell ) const
  -> void
{
  Functor const * func = static_cast< Functor const * >( this );
  func->evaluate( values, coordinates, static_cast< UFCCell const & >( cell ) );
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto
  UserCoefficient< Functor, I, J >::evaluate( size_t            n,
                                              real *            values,
                                              const real *      coordinates,
                                              const ufc::cell & cell ) const
  -> void
{
  Functor const * func = static_cast< Functor const * >( this );

  for ( size_t i = 0; i < n; ++i )
  {
    func->evaluate( &values[i * value_size()],
                    &coordinates[i * cell.geometric_dimension],
                    static_cast< UFCCell const & >( cell ) );
  }
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::eval( real *, const real * ) const
  -> void
{
  error( "UserCoefficient : eval unimplemented" );
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::rank() const
  -> size_t
{
  return ValueSpace< I, J >::rank();
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::dim( size_t i ) const
  -> size_t
{
  return ValueSpace< I, J >::dim( i );
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::value_size() const
  -> size_t
{
  return ValueSpace< I, J >::value_size();
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::interpolate(
  real *                      coefficients,
  const ufc::cell &           ufc_cell,
  const ufc::finite_element & finite_element,
  const Cell &                cell ) const
  -> void
{
  dolfin_assert( coefficients != nullptr );

  if ( ufc_cell_.coordinates.empty() )
  {
    ufc_cell_.init( cell );
  }
  else
  {
    ufc_cell_.update( cell );
  }

  finite_element.evaluate_dofs( coefficients, *this,
                                ufc_cell_.coordinates.data(),
                                ufc_cell_.orientation, ufc_cell );
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::interpolate(
  real *                      coefficients,
  const ufc::cell &           ufc_cell,
  const ufc::finite_element & finite_element,
  const Cell &                cell,
  size_t ) const
  -> void
{
  dolfin_assert( coefficients != nullptr );

  if ( ufc_cell_.coordinates.empty() )
  {
    ufc_cell_.init( cell );
  }
  else
  {
    ufc_cell_.update( cell );
  }

  finite_element.evaluate_dofs( coefficients, *this,
                                ufc_cell_.coordinates.data(),
                                ufc_cell_.orientation, ufc_cell );
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::disp() const
  -> void
{
  section( "UserCoefficient" );
  section( "ValueSpace" );
  message( "rank       : %u", this->rank() );
  message( "value size : %u", this->value_size() );
  end();
  end();
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::sync()
  -> void
{
}

//-----------------------------------------------------------------------------

template < typename Functor, size_t I, size_t J >
inline auto UserCoefficient< Functor, I, J >::sync( Time const & )
  -> void
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_COEFFICENT_H */
