// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Aurélien Larcher 2013-2014.
//
// First added:  2007-04-02
// Last changed: 2014-02-06

#ifndef __DISCRETE_FUNCTION_H
#define __DISCRETE_FUNCTION_H

#include "GenericFunction.h"

#include <dolfin/la/Vector.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>

namespace ufl
{

class FiniteElementBase;

}

namespace dolfin
{

class Mesh;
class Cell;
class Form;
class DofMap;
class FiniteElement;
class SubFunction;
class IntersectionDetector;

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

class DiscreteFunction : public GenericFunction
{
public:

  /// Create discrete function for argument function i of form
  DiscreteFunction(Mesh& mesh, GenericVector& x, Form& form, uint i);

  /// Create discrete function for argument function i of form which owns the vector
  DiscreteFunction(Mesh& mesh, Form& form, uint i);

  /// Create discrete function from given signatures
  DiscreteFunction(Mesh& mesh, GenericVector& x,
                   std::string const& finite_element_signature,
                   std::string const& dof_map_signature);

  /// Create discrete function from given signatures which owns the vector
  DiscreteFunction(Mesh& mesh, std::string const& finite_element_signature,
                   std::string const& dof_map_signature);

  /// Create discrete function from given signatures
  DiscreteFunction(Mesh& mesh, GenericVector& x, std::string const& signature);

  /// Create discrete function from given signatures which owns the vector
  DiscreteFunction(Mesh& mesh, std::string const& signature);

#if ENABLE_UFL
  /// Create discrete function from given UFL Finite Element
  DiscreteFunction(Mesh& mesh, ufl::FiniteElementBase const& finite_element);
#endif

  /// Create discrete function from sub function
  DiscreteFunction(SubFunction& sub_function);

  /// Copy constructor
  DiscreteFunction(DiscreteFunction const& f);

  /// Destructor
  ~DiscreteFunction();

  /// Assign discrete function
  DiscreteFunction const& operator=(const DiscreteFunction& f);

  //--- UFC INTERFACE ---------------------------------------------------------
  /// Evaluate function at given point in cell
  void evaluate(real* values, const real* coordinates,
                const ufc::cell& cell) const;

  //--- GenericFunction -------------------------------------------------------
  /// Return the rank of the value space
  uint rank() const;

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const;

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values(real* values) const;

  /// Interpolate function to finite element space on cell
  void interpolate(real* coefficients, const ufc::cell& cell,
                   const ufc::finite_element& finite_element,
                   const Cell& dolfin_cell) const;

  /// Evaluate function at given point
  void eval(real* values, const real* x) const;

  /// Display basic information
  void disp() const;

  /// Update vector
  void sync_ghosts();

  //---------------------------------------------------------------------------

  /// Return vector
  GenericVector& vector() const;

  /// Return the discrete space
  FiniteElementSpace const& space() const;

  /// Return finite element
  FiniteElement const& finite_element() const;

  /// Return dof map
  DofMap const& dofmap() const;

  /// Return signature
  std::string const signature() const;

  /// Return the number of sub functions
  uint const num_sub_functions() const;

  /// Interpolate values from the given Function
  void interpolate(Function const& other_func);

  /// Get values to cell tabulated block array
  void get_block(real *& values) const;

  /// Set values from cell tabulated block array
  void set_block(real *& values);

  /// Add values from cell tabulated block array
  void add_block(real *& values);

private:

  /// Initialize Vector
  void InitializeVector();

  /// Initialize ghost pattern
  void InitializeGhosts();

  /// Discrete space
  FiniteElementSpace discrete_space_;
  mutable ScratchSpace scratch;

  /// Set to true if local data is owned
  bool const local_vector_;

  /// The vector of dofs
  GenericVector * const X_;

  /// Renumbered dof_map;
  bool renumbered_;
  uint cache_size_;
  uint * indices_;
  real * data_cache_;
  _map<uint, uint> cache_mapping_;

};

}

#endif
