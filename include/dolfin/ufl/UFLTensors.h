// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_TENSORS_H_
#define __UFL_TENSORS_H_

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  ListTensor
 *
 *  @brief  Provides an interface complying with UFL ListTensor.
 *  UFL operator type: Wraps a list of expressions into a tensor
 *  valued expression of one higher rank.
 */

/*
 class ListTensor : public Expression
 {

 public:

 ///
 ListTensor(Expression const& expression);

 ///
 ListTensor(repr_t const & repr);

 ///
 ~ListTensor();

 //--- INTERFACE -------------------------------------------------------------

 ///
 std::vector<Expression const *> const operands() const;

 //      ///
 //      free_indices() const;

 //      ///
 //      index_dimensions() const;

 //      ///
 //      shape() const;

 ///Return the tensor shape of the expression.
 virtual ValueArray const shape() const;

 ///Return a tuple with the free indices (unassigned) of the expression.
 virtual tuple<Index> const free_indices() const;

 ///Return a dict with the free or repeated indices in the expression
 ///as keys and the dimensions of those indices as values.
 virtual dict<IndexBase, type<dolfin::uint> > const index_dimensions();

 /// UFL: Return whether this expression is spatially constant over each cell
 bool const is_cellwise_constant() const;

 //--- INTERFACE inherited from UFLClass -------------------------------------

 /// __repr__
 repr_t const& repr() const;

 /// __str__
 std::string const& str() const;

 ///
 void display() const;

 private:

 ///contains Expression as first and MultiIndex as second
 std::vector<Expression const *> expressions_;
 //      Expression const expression_;
 //      MultiIndex const index_;

 repr_t const repr_;
 std::string const str_;

 };
 */

/**
 *  DOCUMENTATION:
 *
 *  @class  ComponentTensor
 *
 *  @brief  Provides an interface complying with UFL Tensors.
 *  UFL operator type: Maps the free indices of a scalar valued
 *  expression to tensor axes.
 */

class ComponentTensor : public Expression
{

public:

  ///
  ComponentTensor(Expression const& expression, MultiIndex const& multi_index);

  ///
  ComponentTensor(Expression const& expression, tuple<IndexBase> const& indices,
                  dict<IndexBase, type<dolfin::uint> > const& index_dimensions);

  ///
  ComponentTensor(Expression const& expression,
                  tuple<FixedIndex> const& indices);

  ///
  ComponentTensor(repr_t const & repr);

  ///
  ~ComponentTensor();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static ComponentTensor const * create(Object::repr_t const& repr);

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
  bool const is_cellwise_constant() const;

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
  std::vector<Expression const *> const fill_expressions(
      Expression const& e, tuple<IndexBase> const& indices,
      dict<IndexBase, type<dolfin::uint> > const& index_dimensions);
//      std::vector<Expression const *> const fill_expressions(Expression const& e,
//          tuple<FixedIndex> const& indices);

  ///contains Expression as first and MultiIndex as second
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};
} /* namespace ufl */
#endif /* __UFL_TENSORS_H_ */
