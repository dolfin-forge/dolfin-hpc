#ifndef __DOLFIN_FUNCTION_ANALYTIC_H_
#define __DOLFIN_FUNCTION_ANALYTIC_H_

#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

namespace dolfin
{

template < class T >
class Analytic : public GenericFunction
{

public:
  /// Constructor
  Analytic< T >( Mesh & mesh )
    : mesh_( mesh )
    , ufc_cell_( Cell( mesh_, 0 ) )
    , evaluant_()
  {
  }

  /// Constructor
  Analytic< T >( Mesh & mesh, T expr )
    : mesh_( mesh )
    , ufc_cell_( Cell( mesh_, 0 ) )
    , evaluant_( expr )
  {
  }

  /// Destructor
  ~Analytic< T >()
  {
  }

  /// Copy constructor
  Analytic< T >( Analytic< T > const & other )
    : mesh_( const_cast< Mesh & >( other.mesh_ ) )
    , ufc_cell_( other.ufc_cell_ )
    , evaluant_( other.evaluant_ )
  {
  }

  ///
  operator T &()
  {
    return static_cast< T & >( evaluant_ );
  }

  ///
  operator T const &() const
  {
    return static_cast< T const & >( evaluant_ );
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void
    evaluate( real * values, const real * coordinates, const ufc::cell & ) const
  {
    evaluant_.eval( values, coordinates );
  }

  //--- Expression INTERFACE --------------------------------------------------

  ///
  inline void evaluate( size_t            n,
                        real *            values,
                        const real *      coordinates,
                        const ufc::cell & cell ) const
  {
    for ( size_t i = 0; i < n; ++i,
                 values += this->value_size(),
                 coordinates += cell.geometric_dimension )
    {
      evaluant_.eval( values, coordinates );
    }
  }

  /// Evaluate function at given point
  inline void eval( real * values, const real * x ) const
  {
    evaluant_.eval( values, x );
  }

  /// Return the rank of the value space
  inline size_t rank() const
  {
    return evaluant_.rank();
  }

  /// Return the dimension of the value space for axis i
  inline size_t dim( size_t i ) const
  {
    return evaluant_.dim( i );
  }

  /// Value size
  inline size_t value_size() const
  {
    return evaluant_.value_size();
  }

  //--- GenericFunction INTERFACE ---------------------------------------------

  /// Return the mesh
  inline Mesh & mesh() const
  {
    return mesh_;
  }

  /// Interpolate function to vertices of mesh
  inline void interpolate_vertex_values( real * values ) const
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

  /// Interpolate function to finite element space on cell
  inline void interpolate( real *                      coefficients,
                           const ufc::cell &           ufc_cell,
                           const ufc::finite_element & finite_element,
                           const Cell &                cell ) const
  {
    dolfin_assert( coefficients != nullptr );
    ufc_cell_.update( cell );

    finite_element.evaluate_dofs( coefficients, *this,
                                  ufc_cell_.coordinates.data(),
                                  ufc_cell_.orientation, ufc_cell );
  }

  /// Interpolate function to finite element space on facet
  inline void interpolate( real *                      coefficients,
                           const ufc::cell &           ufc_cell,
                           const ufc::finite_element & finite_element,
                           const Cell &                cell,
                           size_t ) const
  {
    dolfin_assert( coefficients != nullptr );
    ufc_cell_.update( cell );

    finite_element.evaluate_dofs( coefficients, *this,
                                  ufc_cell_.coordinates.data(),
                                  ufc_cell_.orientation, ufc_cell );
  }

  /// Display basic information
  inline void disp() const
  {
    section( "Analytic" );
    evaluant_.disp();
    end();
    skip();
  }

  /// Synchronize
  inline void sync()
  {
    // Do nothing
  }

  //---------------------------------------------------------------------------

  /// Delegate time dependency
  Analytic< T > const & operator()( Time const & t ) const
  {
    evaluant_( t );
    return *this;
  }

private:
  /// Evaluant implements time dependency
  void sync( Time const & t )
  {
    evaluant_( t );
  }

  Mesh &          mesh_;
  mutable UFCCell ufc_cell_;
  T               evaluant_;
};

//-----------------------------------------------------------------------------

// Convenience functions to get wrapped evaluant
template < class T >
inline T & evaluant( Analytic< T > & A )
{
  return static_cast< T & >( A );
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_ANALYTIC_H_ */
