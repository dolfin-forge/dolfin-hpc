// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL__TUPLE_H
#define __DOLFIN_UFL__TUPLE_H

#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFL_tuple.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Tuple
 *
 *  @brief  Provides an interface complying with Tuple.
 */

class Tuple : public Expression
{
public:

  ///
  Tuple(tuple<Expression> const& t);

  ///
  Tuple(repr_t const & repr);

  ///
  ~Tuple() override;

  ///
  std::vector<Class const*> const operands(
      std::string const& name) const override;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const override;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  ///
  virtual Tuple const * create(Object::repr_t const& repr) const;

  ///
  std::vector<Expression const *> const operands() const override;

  ///Return the tensor shape of the expression.
  ValueArray const shape() const override;

  ///Return a tuple with the free indices (unassigned) of the expression.
  tuple<Index> const free_indices() const override;

  ///Return a dict with the free or repeated indices in the expression
  ///as keys and the dimensions of those indices as values.
  dict<IndexBase, type<dolfin::uint> > const index_dimensions() const override;

  ///Evaluate the expression tree at the given quadrature_points
  std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const override;

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

//      std::vector<Expression const *> const fill_expressions(std::vector<repr_t> const& reprs);

//      std::vector<Expression const *> const expressions_;
  tuple<Expression> const t_;

  mutable repr_t repr_;
  mutable std::string str_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL__TUPLE_H */
