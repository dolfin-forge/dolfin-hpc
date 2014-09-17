// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Kristian B. Oelgaard, 2007.
// Modified by Martin Sandve Alnes, 2008.
// Modified by Aurélien Larcher 2013. (extension)
//
// First added:  2003-11-28
// Last changed: 2013-06-11

#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/elements/FE.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/common/Variable.h>

#include "SubFunction.h"

#include <ufc.h>

namespace ufl
{

class FiniteElementBase;

}

namespace dolfin
{

class DofMap;
class Expression;
class FiniteElement;
class FiniteElementSpace;
class Form;
class GenericVector;
class Mesh;

/**
 *  @class  Function
 *
 *  @brief  This class represents a function that can be evaluated on a mesh.
 *          The actual representation of the function can vary, but the typical
 *          representation is in terms of a mesh, a vector of degrees of
 *          freedom, a finite element and a dof map that determines the
 *          distribution of degrees of freedom on the mesh.
 *          It is also possible to have user-defined functions, either by
 *          overloading the eval function of this class or by giving a
 *          function (pointer) that returns the value of the function.
 */

class Function : public Variable
{
public:

  /// Function types
  enum Type
  {
    empty = 0, constant, discrete, expression, user
  };

  /// Convert from function type to string
  static std::string type2string(Function::Type type);

  /// Create empty function
  explicit Function(Mesh& mesh);

  /// Create constant scalar function from given value
  explicit Function(Mesh& mesh, real value);

  /// Create constant vector function from given size and value
  explicit Function(Mesh& mesh, uint size, real value);

  /// Create constant vector function from given size and values
  explicit Function(Mesh& mesh, const Array<real>& values);

  /// Create constant tensor function from given shape and values
  explicit Function(Mesh& mesh, const Array<uint>& shape,
                    const Array<real>& values);

  /// Create discrete function for argument function i of form
  explicit Function(Mesh& mesh, GenericVector& x, Form& form, uint i);

  /// Create discrete function for argument function i of form
  explicit Function(Mesh& mesh, Form& form, uint i);

  /// Create discrete function on given discrete space
  Function(GenericVector& x, FiniteElementSpace const& space);

  /// Create discrete function on given discrete space
  Function(FiniteElementSpace const& space);

  /// Create discrete function from UFL Finite Element
  Function(Mesh& mesh, ufl::FiniteElementBase const& finite_element);

  /// Create discrete function on given discrete space [UFC1 compatible]
  Function(Mesh& mesh, std::string const& element, std::string const& dofmap);

  /// Create discrete function from sub function
  explicit Function(SubFunction sub_function);

  /// Assign sub function/slice to discrete function
  Function const& operator=(SubFunction sub_function);

  /// Create expression function
  explicit Function(Mesh& mesh, Expression const& expr);

  /// Assign function
  Function const& operator=(Function& f);

  /// Create function from data file
  explicit Function(const std::string filename);

  /// Copy constructor
  Function(const Function& f);

  /// Destructor
  virtual ~Function();

  //--- DEFERRED INITIALIZATION -----------------------------------------------

  /// Initialize constant function
  void init(Mesh& mesh, real value);

  /// Initialize constant vector-valued function
  void init(Mesh& mesh, uint i, real value);

  /// Initialize discrete function for argument function i of form
  void init(Mesh& mesh, GenericVector& x, Form& form, uint i);

  /// Initialize discrete function for argument function i of form
  void init(Mesh& mesh, Form& form, uint i);

  /// Initialize discrete function on given discrete space
  void init(FiniteElementSpace const& space);

  /// Initialize discrete function on given discrete space
  void init(GenericVector& x, FiniteElementSpace const& space);

  /// Initialize discrete function from UFL Finite Element
  void init(Mesh& mesh, ufl::FiniteElementBase const& finite_element);

  /// Initialize discrete function on given discrete space [UFC1 compatible]
  void init(Mesh& mesh, std::string const& element, std::string const& dofmap);

  /// Initialize expression function
  void init(Mesh& mesh, Expression const& expr);

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  void evaluate(real* values, const real* coordinates,
                const ufc::cell& cell) const;

  //--- COMPOSITION GenericFunction -------------------------------------------

  /// Return the rank of the value space
  virtual uint rank() const;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const;

  /// Return the value size
  virtual uint value_size() const;

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values(real* values) const;

  /// Interpolate function to finite element space on cell
  void interpolate(real* coefficients, const ufc::cell& ufc_cell,
                   const ufc::finite_element& finite_element, Cell& cell,
                   int facet = -1) const;

  /// Evaluate function at given point (overload for user-defined function)
  virtual void eval(real* values, const real* x) const;

  /// Display basic information
  void disp() const;

  /// Synchronize ghosted entries across processes
  void sync_ghosts();

  /// Return the mesh
  Mesh& mesh() const;

  //---------------------------------------------------------------------------

  /// Return the type of function
  Type type() const;

  //--- Wrapper Facade for DiscreteFunction -----------------------------------

  /// Return the vector associated with a DiscreteFunction
  GenericVector& vector() const;

  /// Return the discrete space of a DiscreteFunction
  FiniteElementSpace const& space() const;

  /// Return the finite element space of a DiscreteFunction [TODO: Deprecate]
  FiniteElement const& finite_element() const;

  /// Return the dofmap of a DiscreteFunction [TODO: Deprecate]
  DofMap const& dofmap() const;

  /// Return the signature of a DiscreteFunction
  std::string const signature() const;

  /// Return the number of sub functions of a DiscreteFunction
  uint const num_sub_functions() const;

  /// Interpolate values from the given Function
  void interpolate(const Function& other_func);

  /// Create a cell tabulated block array
  real * create_block() const;

  /// Get values of a DiscreteFunction from cell tabulated block array
  void get_block(real *& values) const;

  /// Set values to a DiscreteFunction from cell tabulated block array
  void set_block(real *& values);

  /// Add values to a DiscreteFunction from cell tabulated block array
  void add_block(real *& values);

  /// Extract sub function/slice from a DiscreteFunction
  SubFunction operator[](uint i);

  //---------------------------------------------------------------------------

protected:

  /// Create user-defined function (evaluation operator must be overloaded)
  explicit Function(Mesh& mesh, uint const& rank, uint const& dim);

  /// Access current cell
  /// (available during assembly for user-defined function)
  Cell const& cell() const;

  /// Access current facet normal
  /// (available during assembly for user-defined function)
  Point normal() const;

  /// Access current facet
  /// (available during assembly for user-defined functions)
  int facet() const;

private:

  /// Default constructor
  Function();

  // Mesh
  Mesh * mesh_;

  // Pointer to current implementation (letter base class)
  GenericFunction* f_;

  // Type of function
  Type type_;

  // Pointer to current cell (if any, otherwise 0)
  mutable Cell* cell_;

  // Current facet (if any, otherwise -1)
  mutable int facet_;

};

//--- INLINES -----------------------------------------------------------------

//--- UFC INTERFACE -----------------------------------------------------------
//-----------------------------------------------------------------------------
inline void Function::evaluate(real* values, const real* coordinates,
                               const ufc::cell& cell) const
{
  f_->evaluate(values, coordinates, cell);
}
//--- COMPOSITION GenericFunction ---------------------------------------------
//-----------------------------------------------------------------------------
inline uint Function::rank() const
{
  return f_->rank();
}
//-----------------------------------------------------------------------------
inline uint Function::dim(unsigned int i) const
{
  return f_->dim(i);
}
//-----------------------------------------------------------------------------
inline uint Function::value_size() const
{
  return f_->value_size();
}
//-----------------------------------------------------------------------------
inline void Function::interpolate_vertex_values(real* values) const
{
  f_->interpolate_vertex_values(values);
}
//-----------------------------------------------------------------------------
inline void Function::interpolate(real* coefficients, const ufc::cell& ufc_cell,
                           const ufc::finite_element& finite_element,
                           Cell& cell, int facet) const
{
  // Make current cell and facet are available to user-defined function
  cell_ = &cell;
  facet_ = facet;

  // Interpolate function
  f_->interpolate(coefficients, ufc_cell, finite_element, cell);

  // Make cell and facet unavailable
  cell_ = 0;
  facet_ = -1;
}
//-----------------------------------------------------------------------------
inline void Function::eval(real* values, const real* x) const
{
  f_->eval(values, x);
}
//-----------------------------------------------------------------------------
inline void Function::sync_ghosts()
{
  f_->sync_ghosts();
}
//-----------------------------------------------------------------------------
inline void Function::disp() const
{
  cout << "Function" << endl;
  cout << "------- " << endl;

  // Begin indentation
  begin("");
  cout << "Type: " << this->type() << " ("
       << Function::type2string(this->type()) << ")" << endl;
  if (f_ != NULL)
  {
    f_->disp();
  }
  // End indentation
  end();
}
//-----------------------------------------------------------------------------
inline Mesh& Function::mesh() const
{
  return f_->mesh();
}
//-----------------------------------------------------------------------------

}

#endif
