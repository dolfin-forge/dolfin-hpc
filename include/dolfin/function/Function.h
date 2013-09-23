// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Kristian B. Oelgaard, 2007.
// Modified by Martin Sandve Alnes, 2008.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2003-11-28
// Last changed: 2013-06-11

#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <dolfin/common/types.h>
#include <dolfin/elements/FE.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/common/Variable.h>

#include "SubFunction.h"

#include <ufc.h>

namespace dolfin
{

  class DofMap;
  class Expression;
  class FiniteElement;
  class Form;
  class GenericFunction;
  class GenericVector;
  class Mesh;

  /// This class represents a function that can be evaluated on a
  /// mesh. The actual representation of the function can vary, but
  /// the typical representation is in terms of a mesh, a vector of
  /// degrees of freedom, a finite element and a dof map that
  /// determines the distribution of degrees of freedom on the mesh.
  ///
  /// It is also possible to have user-defined functions, either by
  /// overloading the eval function of this class or by giving a
  /// function (pointer) that returns the value of the function.

  class Function: public Variable
  {
  public:

    /// Function types
    enum Type
    {
      constant, discrete, empty, expression, ufc, user
    };

    /// Create empty function (read data from file)
    Function();

    /// Create user-defined function (evaluation operator must be overloaded)
    explicit Function(Mesh& mesh);

    /// Create user-defined function from expression
    explicit Function(Mesh& mesh, Expression const& expr);

    /// Create constant scalar function from given value
    Function(Mesh& mesh, real value);

    /// Create constant vector function from given size and value
    Function(Mesh& mesh, uint size, real value);

    /// Create constant vector function from given size and values
    Function(Mesh& mesh, const Array<real>& values);

    /// Create constant tensor function from given shape and values
    Function(Mesh& mesh, const Array<uint>& shape, const Array<real>& values);

    /// Create function from given ufc::function
    Function(Mesh& mesh, const ufc::function& function, uint size);

    /// Create function from given GenericFunction
    //Function(Mesh& mesh, GenericFunction& function);

    /// Create discrete function for argument function i of form
    Function(Mesh& mesh, GenericVector& x, Form& form, uint i = 1);

    /// Create discrete function for argument function i of form
    Function(Mesh& mesh, Form& form, uint i = 1);

    /// Create discrete function from signature
    Function(Mesh& mesh, GenericVector& x,
	     std::string const& finite_element_signature,
	     std::string const& dof_map_signature);

    /// Create discrete function from signature
    Function(Mesh& mesh, std::string const& finite_element_signature,
	     std::string const& dof_map_signature);

    /// Create discrete function from signature
    Function(Mesh& mesh, GenericVector& x,
	     std::string const& finite_element_signature);

    /// Create discrete function from signature
    Function(Mesh& mesh, std::string const& finite_element_signature);

    /// Create discrete function from sub function
    explicit Function(SubFunction& sub_function);

    /// Create function from data file
    explicit Function(const std::string filename);

    /// Copy constructor
    Function(const Function& f);

    /// Destructor
    virtual ~Function();

    /// Create constant function
    void init(Mesh& mesh, real value);

    /// Create expression function
    void init(Mesh& mesh, Expression const& expr);

    /// Create discrete function for argument function i of form
    void init(Mesh& mesh, GenericVector& x, Form& form, uint i = 1);

    /// Create discrete function for argument function i of form
    void init(Mesh& mesh, Form& form, uint i = 1);

    /// Create discrete function from signature
    void init(Mesh& mesh, GenericVector& x,
	      std::string const& finite_element_signature);

    /// Create discrete function from signature
    void init(Mesh& mesh, std::string const& finite_element_signature);

    /// Create discrete function from signature
    void init(Mesh& mesh, std::string const& finite_element_signature,
	      std::string const& dof_map_signature);

    /// Create discrete function from signature
    void init(Mesh& mesh, GenericVector& x,
	      std::string const& finite_element_signature,
	      std::string const& dof_map_signature);

    /// Return the type of function
    Type type() const;

    /// Return the rank of the value space
    virtual uint rank() const;

    /// Return the dimension of the value space for axis i
    virtual uint dim(uint i) const;

    /// Return the mesh
    Mesh& mesh() const;

    //--- !!! BEGIN: Only valid for discrete functions --------------------------

    /// Return the signature of a DiscreteFunction
    std::string signature() const;

    /// Return the vector associated with a DiscreteFunction
    GenericVector& vector() const;

    /// Return the degree of the approximation space of a DiscreteFunction
    uint degree() const;

    /// Return the dofmap of a DiscreteFunction
    DofMap const& dofmap() const;

    /// Return the finite element space of a DiscreteFunction
    FiniteElement const& finite_element() const;

    /// Get values of a DiscreteFunction from array
    void get(real *& values);

    /// Set values to a DiscreteFunction from array
    void set(real *& values);

    /// Return the number of sub functions (only for discrete functions)
    uint numSubFunctions() const;

    /// Extract sub function/slice (only for discrete function)
    SubFunction operator[](uint i);

    //--- !!! END:  Only valid for discrete functions ---------------------------

    /// Assign function
    Function const& operator=(Function& f);

    /// Assign sub function/slice
    Function const& operator=(SubFunction sub_function);

    /// Interpolate function to vertices of mesh
    void interpolate(real* values);

    /// Interpolate values from the given Function
    void interpolate(const Function& other_func);

    /// Interpolate function to finite element space on cell
    void interpolate(real* coefficients, const ufc::cell& ufc_cell,
		     const ufc::finite_element& finite_element, Cell& cell,
		     int facet = -1);

    /// Evaluate function at given point (overload for scalar user-defined function)
    virtual void eval(real* values, const real* x) const;

    /// Evaluate scalar function at given point (overload for scalar user-defined function)
    virtual real eval(const real* x) const;

    void sync_ghosts();

    /// Friends
    friend class LinearPDE;


  protected:

    /// Access current cell (available during assembly for user-defined function)
    Cell const& cell() const;

    /// Access current facet normal (available during assembly for user-defined function)
    Point normal() const;

    /// Access current facet (available during assembly for user-defined functions)
    int facet() const;

  private:

    // Pointer to current implementation (letter base class)
    GenericFunction* f;

    // Type of function
    Type _type;

    // Pointer to current cell (if any, otherwise 0)
    Cell* _cell;

    // Current facet (if any, otherwise -1)
    int _facet;

  };

}

#endif
