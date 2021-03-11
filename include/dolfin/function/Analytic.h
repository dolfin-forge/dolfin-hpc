
#ifndef __DOLFIN_FUNCTION_ANALYTIC_H_
#define __DOLFIN_FUNCTION_ANALYTIC_H_

#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

namespace dolfin
{

template < typename T >
class Analytic : public GenericFunction
{

public:
  /// Constructor
  Analytic( Mesh & mesh );

  /// Constructor
  Analytic( Mesh & mesh, T expr );

  /// Destructor
  ~Analytic();

  /// Copy constructor
  Analytic( Analytic const & other );

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  auto evaluate( real *       values,
                 const real * coordinates,
                 const ufc::cell & ) const -> void override;

  //--- Expression INTERFACE --------------------------------------------------

  ///
  auto evaluate( size_t            n,
                 real *            values,
                 const real *      coordinates,
                 const ufc::cell & cell ) const -> void override;

  /// Evaluate function at given point
  auto eval( real * values, const real * x ) const -> void override;

  /// Return the rank of the value space
  auto rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto dim( size_t i ) const -> size_t override;

  /// Value size
  auto value_size() const -> size_t override;

  //--- GenericFunction INTERFACE ---------------------------------------------

  /// Return the mesh
  auto mesh() const -> Mesh & override;

  /// Interpolate function to vertices of mesh
  auto interpolate_vertex_values( real * values ) const -> void override;

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

  /// Display basic information
  auto disp() const -> void override;

  /// Synchronize
  auto sync() -> void override;

  //---------------------------------------------------------------------------

  /// Delegate time dependency
  auto operator()( Time const & t ) const -> Analytic< T > const &;

private:
  /// Evaluant implements time dependency
  auto sync( Time const & t ) -> void override;

  Mesh & mesh_;
  T      evaluant_;
};

//-----------------------------------------------------------------------------

template < typename T >
inline Analytic< T >::Analytic( Mesh & mesh )
  : mesh_( mesh )
  , evaluant_()
{
}

//-----------------------------------------------------------------------------

template < typename T >
inline Analytic< T >::Analytic( Mesh & mesh, T expr )
  : mesh_( mesh )
  , evaluant_( expr )
{
}

//-----------------------------------------------------------------------------

template < typename T >
inline Analytic< T >::~Analytic()
{
}

//-----------------------------------------------------------------------------

template < typename T >
inline Analytic< T >::Analytic( Analytic< T > const & other )
  : mesh_( const_cast< Mesh & >( other.mesh_ ) )
  , evaluant_( other.evaluant_ )
{
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::evaluate( real *       values,
                                     const real * coordinates,
                                     const ufc::cell & ) const -> void
{
  evaluant_.eval( values, coordinates );
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::evaluate( size_t            n,
                                     real *            values,
                                     const real *      coordinates,
                                     const ufc::cell & cell ) const -> void
{
  for ( size_t i = 0; i < n; ++i, values += this->value_size(),
                                  coordinates += cell.geometric_dimension )
  {
    evaluant_.eval( values, coordinates );
  }
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::eval( real * values, const real * x ) const -> void
{
  evaluant_.eval( values, x );
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::rank() const -> size_t
{
  return evaluant_.rank();
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::dim( size_t i ) const -> size_t
{
  return evaluant_.dim( i );
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::value_size() const -> size_t
{
  return evaluant_.value_size();
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::mesh() const -> Mesh &
{
  return mesh_;
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::interpolate_vertex_values( real * values ) const
  -> void
{
  dolfin_assert( values );

  // Get size of value (number of entries in tensor value)
  size_t const size         = this->value_size();
  real *       local_values = new real[size];
  size_t const num_verts    = mesh_.size( 0 );
  for ( VertexIterator vertex( mesh_ ); !vertex.end(); ++vertex )
  {
    // Evaluate at function at vertex
    evaluant_.eval( local_values, vertex->x() );

    // Copy values to array of vertex values
    for ( size_t i = 0; i < size; ++i )
    {
      values[i * num_verts + vertex->index()] = local_values[i];
    }
  }
  delete[] local_values;
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::interpolate( real *                      coefficients,
                                        UFCCell const &             cell,
                                        ufc::finite_element const & finite_element ) const
  -> void
{
  dolfin_assert( coefficients != nullptr );

  finite_element.evaluate_dofs( coefficients, *this, cell.coordinates.data(),
                                cell.orientation, cell );
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::interpolate( real *                      coefficients,
                                        UFCCell const &             cell,
                                        ufc::finite_element const & finite_element,
                                        size_t ) const -> void
{
  dolfin_assert( coefficients != nullptr );

  finite_element.evaluate_dofs( coefficients, *this, cell.coordinates.data(),
                                cell.orientation, cell );
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::disp() const -> void
{
  section( "Analytic" );
  evaluant_.disp();
  end();
  skip();
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::sync() -> void
{
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::operator()( Time const & t ) const
  -> Analytic< T > const &
{
  evaluant_( t );
  return *this;
}

//-----------------------------------------------------------------------------

template < typename T >
inline auto Analytic< T >::sync( Time const & t ) -> void
{
  evaluant_( t );
}

//-----------------------------------------------------------------------------

// Convenience functions to get wrapped evaluant
template < class T >
inline auto evaluant( Analytic< T > const & A ) -> T const &
{
  return static_cast< T const & >( A );
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_ANALYTIC_H_ */
