// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COEFFICIENT_H
#define __DOLFIN_COEFFICIENT_H

#include <dolfin/common/types.h>
#include <dolfin/evolution/Time.h>

#include <ufc.h>

namespace dolfin
{

class Mesh;
class UFCCell;

/// This class serves as a base class/interface for implementations
/// of specific function representations.

class Coefficient : public ufc::function
{
public:
  /// Constructor
  Coefficient() = default;

  /// Destructor
  ~Coefficient() override = default;

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  void evaluate( real *            values,
                 const real *      coordinates,
                 const ufc::cell & cell ) const override = 0;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given points in cell
  virtual void evaluate( size_t            n,
                         real *            values,
                         const real *      coordinates,
                         const ufc::cell & cell ) const = 0;

  /// Evaluate function at given point
  virtual void eval( real * values, const real * x ) const = 0;

  /// Return the rank of the value space
  virtual auto rank() const -> size_t = 0;

  /// Return the dimension of the value space for axis i
  virtual auto dim( size_t i ) const -> size_t = 0;

  // Return the value size
  virtual auto value_size() const -> size_t = 0;

  /// Interpolate function to finite element space on cell
  virtual void interpolate( real *                      coefficients,
                            UFCCell const &             cell,
                            ufc::finite_element const & finite_element ) const = 0;

  /// Interpolate function to finite element space on facet
  virtual void interpolate( real *                      coefficients,
                            UFCCell const &             cell,
                            ufc::finite_element const & finite_element,
                            size_t                      facet ) const = 0;

  /// Synchronize values
  virtual void sync() = 0;

  /// Display basic information
  virtual void disp() const = 0;

  //---------------------------------------------------------------------------

  /// Delegate time dependence
  auto operator()( Time const & t ) -> Coefficient &
  {
    this->sync( t );
    return *this;
  }

  ///
  template < class T >
  auto compatible( T const & other ) -> bool
  {
    return compatible( *this, other );
  }

  ///
  template < class T >
  static auto compatible( Coefficient const & a, T const & b ) -> bool
  {
    // Check value shape compatibility
    if ( a.rank() != b.rank() )
    {
      warning( "Coefficient: rank mismatch %u != %u", a.rank(), b.rank() );
      return false;
    }
    for ( size_t i = 0; i < a.rank(); ++i )
    {
      if ( a.dim( i ) != b.dim( i ) )
      {
        warning(
          "Coefficient: dim(%u) mismatch %u != %u", i, a.dim( i ), b.dim( i ) );
        return false;
      }
    }
    return true;
  }

private:
  /// Time dependency hook
  virtual void sync( Time const & t ) = 0;
};

//-----------------------------------------------------------------------------

inline auto Coefficient::value_size() const -> size_t
{
  size_t size = 1;
  for ( size_t i = 0; i < this->rank(); ++i )
  {
    size *= this->dim( i );
  }
  return size;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COEFFICIENT_H */
