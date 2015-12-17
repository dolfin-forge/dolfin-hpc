// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_UFL_EXPRESSION_H
#define __DOLFIN_UFL_EXPRESSION_H

#include <dolfin/common/types.h>
#include <dolfin/quadrature/UFCReferenceCell.h>
#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLIntegral.h>
#include <dolfin/ufl/UFLdict.h>
#include <dolfin/ufl/UFL_tuple.h>

using dolfin::uint;

namespace ufl
{
class Argument;
class CoefficientBase;
class FiniteElementBase;
class Index;
class IndexBase;
class UFCReferenceCell;

/**
 *  DOCUMENTATION:
 *
 *  @class  Expression
 *
 *  @brief  Provides an interface complying with UFL Expression.
 */

class Expression : public Class
{
public:

  ///
  ~Expression();

  ///
  static Expression const * create(Object::repr_t const& repr);

  ///
  Expression(std::string const& name);

  ///
  Expression(std::string const& name, repr_t const & repr);

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const = 0;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  virtual std::vector<Expression const *> const operands() const = 0;

  /// UFL: Return the tensor shape of the expression
  virtual ValueArray const shape() const = 0;

  /// UFL: Return the tensor rank of the expression
  uint rank() const;

  /// UFL: Return the cell this expression is defined on
  virtual Cell const cell() const;

  /// UFL: Return the geometric dimension this expression is defined on
  uint geometric_dimension() const;

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  ///Return a tuple with the free indices (unassigned) of the expression.
  virtual tuple<Index> const free_indices() const = 0;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions() const = 0;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const = 0;

  //--- INTERFACE inherited from UFLClass -------------------------------------
  /// __repr__
  virtual repr_t const& repr() const = 0;

  /// __str__
  virtual std::string const& str() const = 0;

  ///
  virtual void display() const = 0;

private:

  bool const is_cellwise_constant_;
};

/*
 class Operator : public Expression
 {
 public:

 Operator (Expression const& expression);
 ~Operator ();

 bool is_cellwise_constant() const;

 private:
 Expression const expression_;
 };
 */
} /* namespace ufl */
#endif /* __DOLFIN_UFL_EXPRESSION_H */
