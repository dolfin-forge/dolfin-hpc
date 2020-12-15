// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_INDEXED_H
#define __DOLFIN_UFL_INDEXED_H

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
 *  @class  Indexed
 *
 *  @brief  Provides an interface complying with UFL Indexed.
 */

class Indexed : public Expression
{

public:

  ///
  Indexed(Expression const& expression, MultiIndex const& multi_index);

  ///
  Indexed(Expression const& expression, tuple<IndexBase> const& indices,
          dict<IndexBase, type<dolfin::uint> > const& index_dimensions);

  ///
//      Indexed(Expression const& expression, tuple<FixedIndex> const& indices);

///
  Indexed(repr_t const & repr);

  ///
  ~Indexed() override;

  ///
  std::vector<Class const*> const operands(
      std::string const& name) const override;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const override;

  //--- INTERFACE -------------------------------------------------------------

  static Indexed const * create(Object::repr_t const& repr);

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
//      std::vector<Expression const *> const fill_expressions(Expression const& e,
//          tuple<FixedIndex> const& indices);

  ///Expression and Index
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};
} /* namespace ufl */
#endif /* __DOLFIN_UFL_INDEXED_H */
