// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_ALGEBRA_H
#define __DOLFIN_UFL_ALGEBRA_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Sum
 *
 *  @brief  Provides an interface complying with UFL Sum.
 */

class Sum : public Expression
{

public:

  ///
  Sum(Expression const& s1, Expression const& s2);

  ///
  Sum(repr_t const& repr);

  ///
  ~Sum();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Sum const * create(Object::repr_t const& repr);

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
  std::vector<Expression const *> const fill_expressions(Expression const& s1,
                                                         Expression const& s2);

  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Product
 *
 *  @brief  Provides an interface complying with UFL Product.
 */

class Product : public Expression
{

public:

  ///
  Product(Expression const& p1, Expression const& p2);

  ///
  Product(repr_t const & repr);

  ///
  ~Product();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Product const * create(Object::repr_t const& repr);

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
  std::vector<Expression const *> const fill_expressions(Expression const& p1,
                                                         Expression const& p2);

  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Division
 *
 *  @brief  Provides an interface complying with UFL Division.
 **/

class Division : public Expression
{

public:

  ///
  Division(Expression const& d1, Expression const& d2);

  ///
  Division(repr_t const & repr);

  ///
  ~Division();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Division const * create(Object::repr_t const& repr);

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
  std::vector<Expression const *> const fill_expressions(Expression const& d1,
                                                         Expression const& d2);

  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Power
 *
 *  @brief  Provides an interface complying with UFL Power.
 */

class Power : public Expression
{

public:

  ///
  Power(Expression const& a, Expression const& b);

  ///
  Power(repr_t const & repr);

  ///
  ~Power();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Power const * create(Object::repr_t const& repr);

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
  std::vector<Expression const *> const fill_expressions(Expression const& a,
                                                         Expression const& b);

  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Abs
 *
 *  @brief  Provides an interface complying with UFL Abs.
 */

class Abs : public Expression
{

public:

  ///
  Abs(Expression const& a);

  ///
  Abs(repr_t const & repr);

  ///
  ~Abs();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Abs const * create(Object::repr_t const& repr);

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
  std::vector<Expression const *> const fill_expressions(Expression const& a);

  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};
} /* namespace ufl */
#endif /* __DOLFIN_UFL_ALGEBRA_H */
