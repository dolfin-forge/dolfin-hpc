// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_DIFFERENTIATION_H
#define __DOLFIN_UFL_DIFFERENTIATION_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLData.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>
#include <dolfin/ufl/UFLTuple.h>
#include <dolfin/ufl/UFLVariable.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  CoefficientDerivative
 *
 *  @brief  Provides an interface complying with UFL CoefficientDerivative.
 *  Derivative of the integrand of a form w.r.t. the degrees of
 *  freedom in a discrete Coefficient.
 */

class CoefficientDerivative : public Expression
{

public:

  ///
  CoefficientDerivative(
      type<dolfin::real> const& integrand, Tuple const& coefficients,
      Tuple const& arguments,
      dict<Data, type<dolfin::uint> > const& coeff_derivatives);

  ///
  CoefficientDerivative(type<dolfin::real> const& integrand,
                        Tuple const& coefficients, Tuple const& arguments,
                        Data const& coeff_derivatives);

  ///
  CoefficientDerivative(repr_t const& repr);

  ///
  ~CoefficientDerivative();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static CoefficientDerivative const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  type<dolfin::real> const integrand_;

  ///two Tuples
  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Tuple const& coefs,
                                                         Tuple const& args);
  std::vector<Expression const *> const expressions_;
  dict<Data, type<dolfin::uint> > const coeff_derivatives_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  SpatialDerivative
 *
 *  @brief  Provides an interface complying with UFL SpatialDerivative.
 */

class SpatialDerivative : public Expression
{

public:

  ///
  SpatialDerivative(Expression const& expression, IndexBase const& index,
                    type<dolfin::uint> const& index_dimension);

  ///
  SpatialDerivative(Expression const& expression, MultiIndex const& index);

  ///
  SpatialDerivative(repr_t const& repr);

  ///
  ~SpatialDerivative();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static SpatialDerivative const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e,
                                                         MultiIndex const& i);
  ///first Expression, second is MultiIndex
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  VariableDerivative
 *
 *  @brief  Provides an interface complying with UFL VariableDerivative.
 */

class VariableDerivative : public Expression
{

public:

  ///
  VariableDerivative(Expression const& expression, Variable const& variable);

  ///
  VariableDerivative(repr_t const& repr);

  ///
  ~VariableDerivative();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static VariableDerivative const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e,
                                                         Variable const& v);
  ///One Expression, one Variable
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  Grad
 *
 *  @brief  Provides an interface complying with UFL Grad.
 */

class Grad : public Expression
{

public:

  ///
  Grad(Expression const& expression);

  ///
  Grad(repr_t const& repr);

  ///
  ~Grad();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Grad const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Div
 *
 *  @brief  Provides an interface complying with UFL Div.
 */

class Div : public Expression
{

public:

  ///
  Div(Expression const& expression);

  ///
  Div(repr_t const& repr);

  ///
  ~Div();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Div const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  NablaGrad
 *
 *  @brief  Provides an interface complying with UFL NablaGrad.
 */

class NablaGrad : public Expression
{

public:

  ///
  NablaGrad(Expression const& expression);

  ///
  NablaGrad(repr_t const& repr);

  ///
  ~NablaGrad();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static NablaGrad const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  NablaDiv
 *
 *  @brief  Provides an interface complying with UFL NablaDiv.
 */

class NablaDiv : public Expression
{

public:

  ///
  NablaDiv(Expression const& expression);

  ///
  NablaDiv(repr_t const& repr);

  ///
  ~NablaDiv();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static NablaDiv const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Curl
 *
 *  @brief  Provides an interface complying with UFL Curl.
 */

class Curl : public Expression
{

public:

  ///
  Curl(Expression const& expression);

  ///
  Curl(repr_t const& repr);

  ///
  ~Curl();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Curl const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  ///Return the tensor shape of the expression.
  virtual ValueArray const shape() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};
} /* namespace ufl */
#endif /* __DOLFIN_UFL_DIFFERENTIATION_H */
