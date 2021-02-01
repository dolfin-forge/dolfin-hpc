
#ifndef __DOLFIN_FUNCTION_REAL_H_
#define __DOLFIN_FUNCTION_REAL_H_

#include <dolfin/fem/Coefficient.h>
#include <dolfin/function/ValueSpace.h>
#include <dolfin/log/log.h>

namespace dolfin
{

template < size_t I = 1, size_t J = 1 >
class Real : public Coefficient
{
public:
  /// Default constructor
  Real()
    : value_()
  {
    std::fill_n( value_, I * J, 0.0 );
  }

  /// Constructor
  explicit Real( real value )
    : value_()
  {
    ( *this ) = value;
  }

  /// Destructor
  ~Real()
  {
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate( real * values, const real *, const ufc::cell & ) const
  {
    std::copy( value_, value_ + I * J, values );
  }

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void
    evaluate( size_t, real * values, const real *, const ufc::cell & ) const
  {
    std::copy( value_, value_ + I * J, values );
  }

  /// Evaluate function at given point
  inline void eval( real * values, const real * ) const
  {
    std::copy( value_, value_ + I * J, values );
  }

  /// Return the rank of the value space
  inline size_t rank() const
  {
    return ValueSpace< I, J >::rank();
  }

  /// Return the dimension of the value space for axis i
  inline size_t dim( size_t i ) const
  {
    return ValueSpace< I, J >::dim( i );
  }

  /// Value size
  inline size_t value_size() const
  {
    return ValueSpace< I, J >::value_size();
  }

  /// Assign constant real number
  inline Real & operator=( real value )
  {
    std::fill_n( value_, I * J, value );
    return *this;
  }

  /// Multiply by constant real number
  inline Real & operator+=( real value )
  {
    for ( size_t i = 0; i < I * J; ++i )
      value_[i] += value;
    return *this;
  }

  /// Add a constant real number
  inline Real & operator-=( real value )
  {
    for ( size_t i = 0; i < I * J; ++i )
      value_[i] -= value;
    return *this;
  }

  /// Substract a constant real number
  inline Real & operator*=( real value )
  {
    for ( size_t i = 0; i < I * J; ++i )
      value_[i] *= value;
    return *this;
  }

  /// Divide by constant real number
  inline Real & operator/=( real value )
  {
    for ( size_t i = 0; i < I * J; ++i )
      value_[i] /= value;
    return *this;
  }

  /// Accessors
  real & operator[]( size_t i )
  {
    return value_[i];
  }
  real const & operator[]( size_t i ) const
  {
    return value_[i];
  }

  ///
  inline Real< I, J > const & operator()( Time const & ) const
  {
    // No-op
    return *this;
  }

  /// Interpolate function to finite element space on cell
  inline void interpolate( real *                      coefficients,
                           const ufc::cell &           ufc_cell,
                           const ufc::finite_element & finite_element,
                           const Cell &                cell ) const
  {
    dolfin_assert( coefficients != nullptr );

    // FIXME this is probably not the smartest way to do it
    size_t const                  gdim     = cell.mesh().geometry_dimension();
    std::vector< size_t > const & vertices = cell.entities( 0 );
    std::vector< double >         coordinates;
    for ( size_t i = 0; i < cell.num_entities( 0 ); ++i )
    {
      double const * coords = cell.mesh().geometry().x( vertices[i] );

      for ( size_t c = 0; c < gdim; ++c )
        coordinates.push_back( coords[c] );
    }

    finite_element.evaluate_dofs(
      coefficients, *this, coordinates.data(), 0, ufc_cell );
  }

  /// Interpolate function to finite element space on facet
  inline void interpolate( real *                      coefficients,
                           const ufc::cell &           ufc_cell,
                           const ufc::finite_element & finite_element,
                           const Cell &                cell,
                           size_t ) const
  {
    dolfin_assert( coefficients != nullptr );

    // FIXME this is probably not the smartest way to do it
    size_t const                  gdim     = cell.mesh().geometry_dimension();
    std::vector< size_t > const & vertices = cell.entities( 0 );
    std::vector< double >         coordinates;
    for ( size_t i = 0; i < cell.num_entities( 0 ); ++i )
    {
      double const * coords = cell.mesh().geometry().x( vertices[i] );

      for ( size_t c = 0; c < gdim; ++c )
        coordinates.push_back( coords[c] );
    }

    finite_element.evaluate_dofs(
      coefficients, *this, coordinates.data(), 0, ufc_cell );
  }

  /// Synchronize
  inline void sync()
  {
    // Do nothing
  }

  /// Display basic information
  inline void disp() const
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

private:
  void sync( Time const & )
  { /* No-op */
  }

  real value_[I * J];
};

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_REAL_H_ */
