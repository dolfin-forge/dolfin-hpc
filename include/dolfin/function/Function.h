// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Aurélien Larcher 2013-2014.
//
// First added:  2007-04-02
// Last changed: 2014-02-06

#ifndef __DOLFIN_FUNCTION_H
#define __DOLFIN_FUNCTION_H

#include <dolfin/function/GenericFunction.h>
#include <dolfin/evolution/TimeDependent.h>

#include <dolfin/common/Array.h>

namespace ufl
{

class FiniteElementSpace;

}

namespace dolfin
{

class Mesh;
class Cell;
class Form;
class DofMap;
class FiniteElement;
class FiniteElementSpace;
class GenericVector;
class ScratchSpace;
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

  /// Create discrete function from given pre-generated UFL Finite Element
  Function(Mesh& mesh, ufl::FiniteElementSpace const& finite_element);

  /// Create discrete function from sub function
  Function(SubFunction const& sub_function);

  /// Copy constructor
  Function(Function const& other);

  /// Destructor
  ~Function();

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
                const ufc::cell& cell) const;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given points in cell
  void evaluate(uint n, real* values, const real* coordinates,
                const ufc::cell& cell) const;

  /// Evaluate function at given point
  void eval(real* values, const real* x) const;

  /// Return the rank of the value space
  uint rank() const;

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const;

  // Return the value size
  uint value_size() const;

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values(real* values) const;

  /// Interpolate function to finite element space on cell
  void interpolate(real* coefficients, const ufc::cell& cell,
                   const ufc::finite_element& finite_element,
                   const Cell& dolfin_cell) const;

  /// Interpolate function to finite element space on facet
  void interpolate(real* coefficients, const ufc::cell& cell,
                   const ufc::finite_element& finite_element,
                   const Cell& dolfin_cell, uint facet) const;

  /// Synchronize values
  void sync();

  /// Display basic information
  void disp() const;

  //---------------------------------------------------------------------------

  /// Return the mesh
  Mesh& mesh() const;

  //---------------------------------------------------------------------------

  /// Return the discrete space
  FiniteElementSpace const& space() const;

  /// Return vector
  GenericVector& vector() const;

  /// Compute interpolation of given function to discrete function
  void interpolate(GenericFunction const& other);

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

  /// AXPY
  Function& axpy(real value, Function const& other);

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
  inline void sync(Time const& t) { TimeDependent::operator ()(t); }

  /// Initialize Vector
  void InitializeVector();

  /// Initialize ghost pattern
  void InitializeGhosts();

  /// Mesh, only allow modification by swap
  Mesh * const mesh_;

  /// Discrete space
  FiniteElementSpace * discrete_space_;
  FiniteElement const * element_;
  DofMap const * dofmap_;
  ScratchSpace * scratch;

  /// Vector of dofs
  GenericVector * X_;

  /// Renumbered dof_map;
  bool renumbered_;
  uint cache_size_;
  uint * indices_;
  real * data_cache_;
  _map<uint, uint> * cache_mapping_;

};

} /* namespace dolfin */

#endif /* DOLFIN_FUNCTION_H */
