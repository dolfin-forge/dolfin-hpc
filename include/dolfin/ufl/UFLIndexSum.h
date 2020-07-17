// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_INDEX_SUM_H
#define __DOLFIN_UFL_INDEX_SUM_H

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>

//#include <dolfin/common/types.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  IndexSum
 *
 *  @brief  Provides an interface complying with UFL IndexSum.
 */

class IndexSum : public Expression
{

public:

  ///
  IndexSum(Expression const& expression, MultiIndex const& multi_index);

  ///
  IndexSum(Expression const& expression, tuple<IndexBase> const& indices,
           dict<IndexBase, type<dolfin::uint> > const& index_dimensions);

  ///
  IndexSum(repr_t const & repr);

  ///
  ~IndexSum() override;

  ///
  std::vector<Class const*> const operands(
      std::string const& name) const override;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const override;

  //--- INTERFACE -------------------------------------------------------------

  static IndexSum const * create(Object::repr_t const& repr);

  ///
  MultiIndex const * index() const;

  ///
  dolfin::uint dimension() const;

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

  /// UFL: Return whether this expression is spatially constant over each cell
  bool is_cellwise_constant() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e,
                                                         MultiIndex const& i);
  std::vector<Expression const *> const fill_expressions(
      Expression const& e, tuple<IndexBase> const& indices,
      dict<IndexBase, type<dolfin::uint> > const& index_dimensions);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};
} /* namespace ufl */
#endif /* __DOLFIN_UFL_INDEX_SUM_H */
