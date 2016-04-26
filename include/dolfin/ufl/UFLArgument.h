// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_UFL_ARGUMENT_H
#define __DOLFIN_UFL_ARGUMENT_H

#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLFiniteElement.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Argument
 *
 *  @brief  Provides an interface complying with UFL Argument.
 */

class Argument : public Expression
{
public:

  ///
  Argument(FiniteElementSpace const& finite_element, dolfin::uint const& count);

  ///
  Argument(repr_t const& repr);

  ///
  ~Argument();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Argument const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  /// Return a reference to the FiniteElementSpace of this Argument
  FiniteElementSpace const& element() const;

  /// Return the count of this Argument
  type<dolfin::uint> const& count() const;

  /// Return a reference to the cell of the FiniteElementSpace of this Argument
  Cell const cell() const;

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  bool is_cellwise_constant() const;

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

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

protected:

  FiniteElementSpace const& finite_element_;
  type<dolfin::uint> const count_;

//      std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_ARGUMENT_H */
