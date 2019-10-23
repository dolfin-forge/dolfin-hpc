// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_COEFFICIENT_H
#define __DOLFIN_UFL_COEFFICIENT_H

#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLFiniteElement.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  CoefficientBase
 *
 *  @brief  Base class for interface complying with UFL Coefficient.
 */

class CoefficientBase : public Expression
{

public:

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const = 0;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  static CoefficientBase const * create(Object::repr_t const& repr);

  ///
  std::vector<Expression const *> const operands() const;

  /// Return a reference to the FiniteElementSpace of this Coefficient
  FiniteElementSpace const& element() const;

  /// Return a reference to the cell of the FiniteElementSpace of this Coefficient
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
  repr_t const& repr() const = 0;

  /// __str__
  std::string const& str() const = 0;

  ///
  void display() const = 0;

protected:

  ///
  CoefficientBase(std::string const& name,
                  FiniteElementSpace const& finite_element,
                  dolfin::uint const& count);

  ///
  CoefficientBase(std::string const& name, Cell const& cell,
                  dolfin::uint const& count);

  ///
  CoefficientBase(std::string const& name, Cell const& cell,
                  dolfin::uint const& dim, dolfin::uint const& count);

  ///
  CoefficientBase(std::string const& name, Cell const& cell,
                  ValueArray const& shape,
                  std::map<dolfin::uint, dolfin::uint> const& symmetry,
                  dolfin::uint const& count);

  ///
  CoefficientBase(std::string const& name, repr_t const& repr);

  ///
  ~CoefficientBase();

  FiniteElementSpace const* finite_element_;
  type<dolfin::uint> const count_;

//      std::vector<Expression const *> const expressions_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  Coefficient
 *
 *  @brief  Provides an interface complying with UFL Coefficient.
 */

class Coefficient : public CoefficientBase
{
public:

  ///
  Coefficient(FiniteElementSpace const& finite_element,
              dolfin::uint const& count);

  ///
  Coefficient(repr_t const& repr);

  ///
  ~Coefficient();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  static Coefficient const * create(Object::repr_t const& repr);

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

protected:

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  Constant
 *
 *  @brief  Provides an interface complying with UFL Constant.
 */

class Constant : public CoefficientBase
{
public:

  ///
  Constant(Cell const& cell, dolfin::uint const& count);

  ///
  Constant(repr_t const& repr);

  ///
  ~Constant();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const*> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static Constant const * create(Object::repr_t const& repr);

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

protected:

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  VectorConstant
 *
 *  @brief  Provides an interface complying with UFL VectorConstant.
 */

class VectorConstant : public CoefficientBase
{
public:

  ///
  VectorConstant(Cell const& cell, dolfin::uint const& dim,
                 dolfin::uint const& count);

  ///
  VectorConstant(repr_t const& repr);

  ///
  ~VectorConstant();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const*> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static VectorConstant const * create(Object::repr_t const& repr);

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

protected:

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  TensorConstant
 *
 *  @brief  Provides an interface complying with UFL TensorConstant.
 */

class TensorConstant : public CoefficientBase
{
public:

  ///
  TensorConstant(Cell const& cell, ValueArray const& shape,
                 std::map<dolfin::uint, dolfin::uint> const& symmetry,
                 dolfin::uint const& count);

  ///
  TensorConstant(repr_t const& repr);

  ///
  ~TensorConstant();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const*> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static TensorConstant const * create(Object::repr_t const& repr);

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

protected:

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_COEFFICIENT_H */
