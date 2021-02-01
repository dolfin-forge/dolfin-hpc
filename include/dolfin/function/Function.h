// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTION_H
#define __DOLFIN_FUNCTION_H

#include <dolfin/evolution/TimeDependent.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/FunctionDecomposition.h>
#include <dolfin/function/FunctionInterpolation.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/la/GenericVector.h>

#include <vector>

namespace ufl
{

class FiniteElementSpace;

}

namespace dolfin
{

class Cell;
class Expression;
class Form;
class Mesh;
class SubFunction;

/**
 *  @class  DiscreteFunction
 *
 *  @brief  This class implements the functionality for discrete functions.
 *          A discrete function is defined in terms of a mesh, a vector of
 *          degrees of freedom, a finite element and a dof map. The finite
 *          element determines how the function is defined locally on each
 *          cell of the mesh in terms of the local degrees of freedom, and
 *          the dof map determines how the degrees of freedom are
 *          distributed on the mesh.
 */

class Function : public GenericFunction, public TimeDependent
{
public:
  /// Create discrete function for argument function i of form
  /// The discrete space is defined on the i-th coefficient mesh.
  Function( Form & form, size_t i );

  /// Create discrete function from given  discrete space
  Function( FiniteElementSpace const & space );

  /// Create discrete function from sub function
  Function( SubFunction const & sub_function );

  /// Copy constructor
  Function( Function const & other );

  /// Destructor
  ~Function() override;

  //--- DEFERRED INITIALIZATION -----------------------------------------------
  // NOTE: Beware, camembert !

  /// Create an empty discrete function without any mesh assignment
  Function();

  /// Create an empty discrete function
  Function( Mesh & mesh );

  /// Return whether the function is empty
  auto empty() const -> bool;

  /// Initialize discrete function for argument function i of form
  void init( Form & form, size_t i );

  /// Initialize discrete function on given discrete space
  void init( FiniteElementSpace const & space );

  /// Clear attributes, function becomes empty
  void clear();

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  void evaluate( real *            values,
                 const real *      coordinates,
                 const ufc::cell & cell ) const override;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given points in cell
  void evaluate( size_t            n,
                 real *            values,
                 const real *      coordinates,
                 const ufc::cell & cell ) const override;

  /// Evaluate function at given point
  void eval( real * values, const real * x ) const override;

  /// Return the rank of the value space
  auto rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto dim( size_t i ) const -> size_t override;

  // Return the value size
  auto value_size() const -> size_t override;

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values( real * values ) const override;

  /// Interpolate function to finite element space on cell
  void interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell &                dolfin_cell ) const override;

  /// Interpolate function to finite element space on facet
  void interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell &                dolfin_cell,
                    size_t                      facet ) const override;

  /// Synchronize values
  void sync() override;

  /// Display basic information
  void disp() const override;

  //---------------------------------------------------------------------------

  /// Return the mesh
  auto mesh() const -> Mesh & override;

  //---------------------------------------------------------------------------

  /// Return the discrete space
  auto space() const -> FiniteElementSpace const &;

  /// Return vector
  auto vector() const -> GenericVector &;

  /// Interpolate from an expression
  void operator<<( Expression const & other );

  /// Interpolate from a coefficient
  void operator<<( Coefficient const & other );

  /// Interpolate from a generic function
  void operator<<( GenericFunction const & other );

  /// Compute function decomposition into scalar component functions
  auto decompose() -> std::vector< Function * >;

  /// Return the number of sub functions i.e number of subspaces
  auto num_sub_functions() const -> size_t;

  /// Get the size of tabulated block array
  auto block_size() const -> uidx;

  /// Create a cell tabulated block array
  auto create_block() const -> real *;

  /// Get values to cell tabulated block array
  void get_block( real *& values ) const;

  /// Set values from cell tabulated block array
  void set_block( real *& values );

  /// Add values from cell tabulated block array
  void add_block( real *& values );

  /// Swap functions
  auto swap( Function & other ) -> Function &;

  //--- OPERATORS -------------------------------------------------------------
  // NOTE: A sequence of operations on Functions must be finalized with a sync.

  /// Assignment, fill the function if empty, error if space mismatch
  auto operator=( Function const & other ) -> Function &;

  /// Addition, error if empty or space mismatch
  auto operator+=( Function const & other ) -> Function &;

  /// Subtraction, error if empty or space mismatch
  auto operator-=( Function const & other ) -> Function &;

  /// Component-wise multiplication, error if empty or space mismatch
  auto operator*=( Function const & other ) -> Function &;

  /// Add multiple of given vector (y=a*x+y)
  auto axpy( real a, const Function & x ) -> Function &;

  /// Add multiple of given vector (y=a*x+b*y)
  auto axpby( real a, const Function & x, real b ) -> Function &;

  /// Add multiple of given vector (w=a*x+y)
  auto waxpy( real a, const Function & x, const Function & y ) -> Function &;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  auto
    axpbypcz( real a, const Function & x, real b, const Function & y, real c )
      -> Function &;

  /// Assignment to a real value, error if empty
  auto operator=( real value ) -> Function &;

  /// Addition to a real number, error if empty [Not implemented]
  auto operator+=( real value ) -> Function &;

  /// Subtraction of a real number, error if empty [Not implemented]
  auto operator-=( real value ) -> Function &;

  /// Multiplication with a real number, error if empty
  auto operator*=( real value ) -> Function &;

  /// Division by a real number, error if empty
  auto operator/=( real value ) -> Function &;

  /// Zero the vector
  auto zero() -> Function &;

  /// Return the minimum value
  auto min() const -> real;

  /// Return the maximum value
  auto max() const -> real;

  //---------------------------------------------------------------------------

  /// Time dependency
  auto operator()( Time const & t ) -> Function &
  {
    TimeDependent::operator()( t );
    return *this;
  }

private:
  /// Time synchronization hook
  inline void sync( Time const & t ) override
  {
    TimeDependent::operator()( t );
  }

  /// Initialize Vector
  void InitializeVector();

  /// Initialize ghost pattern
  void InitializeGhosts();

private:
  /// Mesh, only allow modification by swap or assignment
  Mesh * const mesh_ { nullptr };

  /// FiniteElement space
  FiniteElementSpace * fe_space_ { nullptr };

  // temporary storage
  ScratchSpace * scratch_ { nullptr };

  /// Vector of dofs
  GenericVector * X_ { nullptr };

  /// Renumbered dof_map;
  bool renumbered_ { false };

  /// function cache
  std::vector< size_t > indices_;
  std::vector< real >   data_cache_;

  /// cached mappings
  _map< size_t, size_t >         cache_mapping_;
  _ordered_map< size_t, size_t > ghost_mapping_;
};

//-----------------------------------------------------------------------------

inline auto Function::empty() const -> bool
{
  return ( fe_space_ == nullptr );
}

//--- UFC INTERFACE -----------------------------------------------------------

inline void Function::evaluate( real *            values,
                                const real *      x,
                                const ufc::cell & cell ) const
{
  evaluate( 1, values, x, cell );
}

//--- GenericFunction ---------------------------------------------------------

inline auto Function::mesh() const -> Mesh &
{
  return ( *mesh_ );
}

//-----------------------------------------------------------------------------

inline auto Function::rank() const -> size_t
{
  dolfin_assert( fe_space_ != nullptr );
  return fe_space_->element().value_rank();
}

//-----------------------------------------------------------------------------

inline auto Function::dim( size_t i ) const -> size_t
{
  dolfin_assert( fe_space_ != nullptr );
  return fe_space_->element().value_dimension( i );
}

//-----------------------------------------------------------------------------

inline auto Function::value_size() const -> size_t
{
  dolfin_assert( scratch_ != nullptr );
  return scratch_->size;
}
//-----------------------------------------------------------------------------

inline void Function::interpolate( real *                      coefficients,
                                   const ufc::cell &           cell,
                                   const ufc::finite_element & finite_element,
                                   const Cell &                dolfin_cell,
                                   size_t ) const
{
  interpolate( coefficients, cell, finite_element, dolfin_cell );
}

//-----------------------------------------------------------------------------

inline auto Function::vector() const -> GenericVector &
{
  dolfin_assert( X_ );
  return *X_;
}

//-----------------------------------------------------------------------------

inline auto Function::space() const -> FiniteElementSpace const &
{
  dolfin_assert( fe_space_ );
  return *fe_space_;
}

//-----------------------------------------------------------------------------

inline void Function::operator<<( Expression const & other )
{
  FunctionInterpolation::compute( other, *this );
}

//-----------------------------------------------------------------------------

inline void Function::operator<<( Coefficient const & other )
{
  FunctionInterpolation::compute( other, *this );
}

//-----------------------------------------------------------------------------

inline void Function::operator<<( GenericFunction const & other )
{
  FunctionInterpolation::compute( other, *this );
}

//-----------------------------------------------------------------------------

inline auto Function::decompose() -> std::vector< Function * >
{
  return FunctionDecomposition::compute( *this );
}

//-----------------------------------------------------------------------------

inline auto Function::num_sub_functions() const -> size_t
{
  dolfin_assert( fe_space_ != nullptr );
  return fe_space_->element().num_sub_elements();
}

//-----------------------------------------------------------------------------

inline auto Function::block_size() const -> uidx
{
  dolfin_assert( fe_space_ != nullptr );
  return fe_space_->dofmap().dofsmapping_size();
}

//-----------------------------------------------------------------------------

inline auto Function::create_block() const -> real *
{
  dolfin_assert( fe_space_ != nullptr );
  return new real[fe_space_->dofmap().dofsmapping_size()];
}

//-----------------------------------------------------------------------------

inline void Function::get_block( real *& values ) const
{
  dolfin_assert( X_ != nullptr );
  dolfin_assert( fe_space_ != nullptr );
  if ( !values )
  {
    values = new real[fe_space_->dofmap().dofsmapping_size()];
  }
  X_->apply();
  X_->get( values,
           fe_space_->dofmap().dofsmapping_size(),
           fe_space_->dofmap().dofsmapping() );
}

//-----------------------------------------------------------------------------

inline void Function::set_block( real *& values )
{
  dolfin_assert( X_ != nullptr );
  dolfin_assert( fe_space_ != nullptr );
  X_->set( values,
           fe_space_->dofmap().dofsmapping_size(),
           fe_space_->dofmap().dofsmapping() );
  sync();
}

//-----------------------------------------------------------------------------

inline void Function::add_block( real *& values )
{
  dolfin_assert( X_ != nullptr );
  dolfin_assert( fe_space_ != nullptr );
  X_->add( values,
           fe_space_->dofmap().dofsmapping_size(),
           fe_space_->dofmap().dofsmapping() );
  sync();
}

//-----------------------------------------------------------------------------

inline auto Function::operator+=( Function const & other ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() += other.vector();
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator-=( Function const & other ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() -= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator*=( Function const & other ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() *= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::axpy( real a, const Function & x ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( this->space() == x.space() );
  this->vector().axpy( a, x.vector() );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::axpby( real a, const Function & x, real b ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( this->space() == x.space() );
  this->vector().axpby( a, x.vector(), b );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::waxpy( real a, const Function & x, const Function & y )
  -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( !y.empty() );
  dolfin_assert( this->space() == x.space() );
  dolfin_assert( this->space() == y.space() );
  this->vector().waxpy( a, x.vector(), y.vector() );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::axpbypcz( real             a,
                                const Function & x,
                                real             b,
                                const Function & y,
                                real             c ) -> Function &
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( !y.empty() );
  dolfin_assert( this->space() == x.space() );
  dolfin_assert( this->space() == y.space() );
  this->vector().axpbypcz( a, x.vector(), b, y.vector(), c );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator=( real value ) -> Function &
{
  dolfin_assert( !this->empty() );
  this->vector() = value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator+=( real ) -> Function &
{
  dolfin_assert( !this->empty() );
  error( "Not implemented" );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator-=( real ) -> Function &
{
  dolfin_assert( !this->empty() );
  error( "Not implemented" );
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator*=( real value ) -> Function &
{
  dolfin_assert( !this->empty() );
  this->vector() *= value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::operator/=( real value ) -> Function &
{
  dolfin_assert( !this->empty() );
  this->vector() /= value;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::zero() -> Function &
{
  dolfin_assert( !this->empty() );
  this->vector().zero();
  return *this;
}

//-----------------------------------------------------------------------------

inline auto Function::min() const -> real
{
  dolfin_assert( !this->empty() );
  return this->vector().min();
}

//-----------------------------------------------------------------------------

inline auto Function::max() const -> real
{
  dolfin_assert( !this->empty() );
  return this->vector().max();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* DOLFIN_FUNCTION_H */
