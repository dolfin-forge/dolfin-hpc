// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTION_H
#define __DOLFIN_FUNCTION_H

#include <dolfin/common/Array.h>
#include <dolfin/evolution/TimeDependent.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/FunctionDecomposition.h>
#include <dolfin/function/FunctionInterpolation.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/la/GenericVector.h>

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
  Function(Form& form, uint i);

  /// Create discrete function from given  discrete space
  Function(FiniteElementSpace const& space);

  /// Create discrete function from sub function
  Function(SubFunction const& sub_function);

  /// Copy constructor
  Function(Function const& other);

  /// Destructor
  ~Function() override;

  //--- DEFERRED INITIALIZATION -----------------------------------------------
  // NOTE: Beware, camembert !

  /// Create an empty discrete function without any mesh assignment
  Function();

  /// Create an empty discrete function
  Function(Mesh& mesh);

  /// Return whether the function is empty
  bool empty() const;

  /// Initialize discrete function for argument function i of form
  void init(Form& form, uint i);

  /// Initialize discrete function on given discrete space
  void init(FiniteElementSpace const& space);

  /// Clear attributes, function becomes empty
  void clear();

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  void evaluate(real* values, const real* coordinates,
                const ufc::cell& cell) const override;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given points in cell
  void evaluate(uint n, real* values, const real* coordinates,
                const ufc::cell& cell) const override;

  /// Evaluate function at given point
  void eval(real* values, const real* x) const override;

  /// Return the rank of the value space
  uint rank() const override;

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const override;

  // Return the value size
  uint value_size() const override;

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values(real* values) const override;

  /// Interpolate function to finite element space on cell
  void interpolate(real* coefficients, const ufc::cell& cell,
                   const ufc::finite_element& finite_element,
                   const Cell& dolfin_cell) const override;

  /// Interpolate function to finite element space on facet
  void interpolate(real* coefficients, const ufc::cell& cell,
                   const ufc::finite_element& finite_element,
                   const Cell& dolfin_cell, uint facet) const override;

  /// Synchronize values
  void sync() override;

  /// Display basic information
  void disp() const override;

  //---------------------------------------------------------------------------

  /// Return the mesh
  Mesh& mesh() const override;

  //---------------------------------------------------------------------------

  /// Return the discrete space
  FiniteElementSpace const& space() const;

  /// Return vector
  GenericVector& vector() const;

  /// Interpolate from an expression
  void operator<<(Expression const& other);

  /// Interpolate from a coefficient
  void operator<<(Coefficient const& other);

  /// Interpolate from a generic function
  void operator<<(GenericFunction const& other);

  /// Compute function decomposition into scalar component functions
  Array<Function *> decompose();

  /// Return the number of sub functions i.e number of subspaces
  uint num_sub_functions() const;

  /// Get the size of tabulated block array
  uidx block_size() const;

  /// Create a cell tabulated block array
  real * create_block() const;

  /// Get values to cell tabulated block array
  void get_block(real *& values) const;

  /// Set values from cell tabulated block array
  void set_block(real *& values);

  /// Add values from cell tabulated block array
  void add_block(real *& values);

  /// Swap functions
  Function& swap(Function& other);

  //--- OPERATORS -------------------------------------------------------------
  // NOTE: A sequence of operations on Functions must be finalized with a sync.

  /// Assignment, fill the function if empty, error if space mismatch
  Function& operator=(Function const& other);

  /// Addition, error if empty or space mismatch
  Function& operator+=(Function const& other);

  /// Subtraction, error if empty or space mismatch
  Function& operator-=(Function const& other);

  /// Component-wise multiplication, error if empty or space mismatch
  Function& operator*=(Function const& other);

  /// Add multiple of given vector (y=a*x+y)
  Function& axpy( real a, const Function & x );

  /// Add multiple of given vector (y=a*x+b*y)
  Function& axpby( real a, const Function & x,
                   real b );

  /// Add multiple of given vector (w=a*x+y)
  Function& waxpy( real a, const Function & x,
                           const Function & y );

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  Function& axpbypcz( real a, const Function & x,
                      real b, const Function & y,
                      real c );

  /// Assignment to a real value, error if empty
  Function& operator=(real value);

  /// Addition to a real number, error if empty [Not implemented]
  Function& operator+=(real value);

  /// Subtraction of a real number, error if empty [Not implemented]
  Function& operator-=(real value);

  /// Multiplication with a real number, error if empty
  Function& operator*=(real value);

  /// Division by a real number, error if empty
  Function& operator/=(real value);

  /// Zero the vector
  Function& zero();

  /// Return the minimum value
  real min() const;

  /// Return the maximum value
  real max() const;

  //---------------------------------------------------------------------------

  /// Time dependency
  Function& operator()(Time const& t)
  {
    TimeDependent::operator ()(t);
    return *this;
  }

private:

  /// Time synchronization hook
  inline void sync(Time const& t) override { TimeDependent::operator ()(t); }

  /// Initialize Vector
  void InitializeVector();

  /// Initialize ghost pattern
  void InitializeGhosts();

  /// Mesh, only allow modification by swap or assignment
  Mesh * const mesh_{nullptr};

  /// Discrete space
  FiniteElementSpace * discrete_space_{nullptr};
  // FIXME element_ and dofmap_ are already part of the discrete_space_
  FiniteElement const * element_{nullptr};
  DofMap const * dofmap_{nullptr};
  ScratchSpace * scratch{nullptr};

  /// Vector of dofs
  GenericVector * X_{nullptr};

  /// Renumbered dof_map;
  bool renumbered_{false};
  uint cache_size_{0};
  uint * indices_{nullptr};
  real * data_cache_{nullptr};
  _map<uint, uint> * cache_mapping_{nullptr};

};

//-----------------------------------------------------------------------------
inline bool Function::empty() const
{
  return ( discrete_space_ == nullptr );
}

//--- UFC INTERFACE -----------------------------------------------------------
inline void Function::evaluate( real *            values,
                         const real *      x,
                         const ufc::cell & cell ) const
{
  evaluate( 1, values, x, cell );
}

//--- GenericFunction ---------------------------------------------------------
inline Mesh & Function::mesh() const
{
  return ( *mesh_ );
}

//-----------------------------------------------------------------------------
inline uint Function::rank() const
{
  dolfin_assert( element_ );
  return element_->value_rank();
}

//-----------------------------------------------------------------------------
inline uint Function::dim( uint i ) const
{
  dolfin_assert( element_ );
  return element_->value_dimension( i );
}

//-----------------------------------------------------------------------------
inline uint Function::value_size() const
{
  dolfin_assert( scratch );
  return scratch->size;
}
//-----------------------------------------------------------------------------
inline void Function::interpolate( real *                      coefficients,
                            const ufc::cell &           cell,
                            const ufc::finite_element & finite_element,
                            const Cell &                dolfin_cell,
                            uint ) const
{
  interpolate( coefficients, cell, finite_element, dolfin_cell );
}

//-----------------------------------------------------------------------------
inline GenericVector & Function::vector() const
{
  dolfin_assert( X_ );
  return *X_;
}

//-----------------------------------------------------------------------------
inline FiniteElementSpace const & Function::space() const
{
  dolfin_assert( discrete_space_ );
  return *discrete_space_;
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
inline Array< Function * > Function::decompose()
{
  return FunctionDecomposition::compute( *this );
}

//-----------------------------------------------------------------------------
inline uint Function::num_sub_functions() const
{
  dolfin_assert( element_ );
  return element_->num_sub_elements();
}

//-----------------------------------------------------------------------------
inline uidx Function::block_size() const
{
  dolfin_assert( dofmap_ );
  return dofmap_->dofsmapping_size();
}

//-----------------------------------------------------------------------------
inline real * Function::create_block() const
{
  dolfin_assert( dofmap_ );
  return new real[dofmap_->dofsmapping_size()];
}

//-----------------------------------------------------------------------------
inline void Function::get_block( real *& values ) const
{
  dolfin_assert( X_ );
  dolfin_assert( dofmap_ );
  if ( !values )
  {
    values = new real[dofmap_->dofsmapping_size()];
  }
  X_->apply();
  X_->get( values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping() );
}

//-----------------------------------------------------------------------------
inline void Function::set_block( real *& values )
{
  dolfin_assert( X_ );
  dolfin_assert( dofmap_ );
  X_->set( values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping() );
  sync();
}

//-----------------------------------------------------------------------------
inline void Function::add_block( real *& values )
{
  dolfin_assert( X_ );
  dolfin_assert( dofmap_ );
  X_->add( values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping() );
  sync();
}

//-----------------------------------------------------------------------------
inline Function & Function::operator+=( Function const & other )
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() += other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator-=( Function const & other )
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() -= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator*=( Function const & other )
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !other.empty() );
  dolfin_assert( this->space() == other.space() );
  this->vector() *= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::axpy( real a, const Function & x )
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( this->space() == x.space() );
  this->vector().axpy( a, x.vector() );
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::axpby( real a, const Function & x,
                                   real b )
{
  dolfin_assert( !this->empty() );
  dolfin_assert( !x.empty() );
  dolfin_assert( this->space() == x.space() );
  this->vector().axpby( a, x.vector(), b );
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::waxpy( real a, const Function & x,
                                           const Function & y )
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
inline Function & Function::axpbypcz( real a, const Function & x,
                                      real b, const Function & y,
                                      real c )
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
inline Function & Function::operator=( real value )
{
  dolfin_assert( !this->empty() );
  this->vector() = value;
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator+=( real )
{
  dolfin_assert( !this->empty() );
  error( "Not implemented" );
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator-=( real )
{
  dolfin_assert( !this->empty() );
  error( "Not implemented" );
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator*=( real value )
{
  dolfin_assert( !this->empty() );
  this->vector() *= value;
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::operator/=( real value )
{
  dolfin_assert( !this->empty() );
  this->vector() /= value;
  return *this;
}

//-----------------------------------------------------------------------------
inline Function & Function::zero()
{
  dolfin_assert( !this->empty() );
  this->vector().zero();
  return *this;
}

//-----------------------------------------------------------------------------
inline real Function::min() const
{
  dolfin_assert( !this->empty() );
  return this->vector().min();
}

//-----------------------------------------------------------------------------
inline real Function::max() const
{
  dolfin_assert( !this->empty() );
  return this->vector().max();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* DOLFIN_FUNCTION_H */
