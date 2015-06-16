// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_INDEX_H_
#define __UFL_INDEX_H_

#include <dolfin/ufl/UFLdict.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLtuple.h>
#include <dolfin/ufl/UFLtype.h>

#include <dolfin/common/types.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  IndexBase
 *
 *  @brief  Provides an interface complying with UFL IndexBase.
 */

class IndexBase : public Expression
{

public:

//      ///
//      IndexBase(dolfin::uint const& count);

///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const = 0;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static IndexBase const * create(Object::repr_t const& repr);

  ///
  type<dolfin::uint> const & count() const;

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
  virtual repr_t const& repr() const = 0;

  /// __str__
  virtual std::string const& str() const = 0;

  ///
  virtual void display() const = 0;

protected:

  ///
  IndexBase(std::string const& name, dolfin::uint const& count);

  ///
  IndexBase(std::string const& name, IndexBase const& index);

  ///
  IndexBase(std::string const & name, repr_t const & repr);

  ///
  virtual ~IndexBase();

//    private:

  type<dolfin::uint> const count_;
  std::vector<Expression const *> const expressions_;

private:
  std::string const name_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Index
 *
 *  @brief  Provides an interface complying with UFL Index.
 *  UFL value: An index with no value assigned.
 *  Used to represent free indices in Einstein indexing notation.
 */

class Index : public IndexBase
{

public:

  ///
  Index(type<dolfin::uint> const& count);

  ///
  Index(repr_t const & repr);

  ///
  ~Index();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static Index const * create(Object::repr_t const& repr);

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  FixedIndex
 *
 *  @brief  Provides an interface complying with UFL FixedIndex.
 *  UFL value: An index with a specific value assigned.
 */

class FixedIndex : public IndexBase
{

public:

  ///
  FixedIndex(dolfin::uint const& value);

  ///
  FixedIndex(repr_t const & repr);

  ///
  ~FixedIndex();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static FixedIndex const * create(Object::repr_t const& repr);

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  repr_t const repr_;
  std::string const str_;

};

/**
 *  DOCUMENTATION:
 *
 *  @class  MultiIndex
 *
 *  @brief  Provides an interface complying with UFL MultiIndex.
 *  Represents a sequence of indices, either fixed or free.
 */

class MultiIndex : public Expression
{

public:

  ///
//      MultiIndex(tuple<type<dolfin::uint> > const& indices,
//          dict<type<dolfin::uint>, type<dolfin::uint> > const& index_dimensions);

  ///
  MultiIndex(tuple<IndexBase> const& indices,
             dict<IndexBase, type<dolfin::uint> > const& index_dimensions);

//      ///
//      MultiIndex(tuple<FixedIndex> const& indices);

  ///
  MultiIndex(repr_t const & repr);

  ///
  ~MultiIndex();

  ///
  virtual std::vector<Class const*> const operands(
      std::string const& name) const;

  ///
  virtual std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  //--- INTERFACE -------------------------------------------------------------

  ///
  static MultiIndex const * create(Object::repr_t const& repr);

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

  tuple<IndexBase> const fill_indices(repr_t const& repr);
//      tuple<FixedIndex> const fill_fixed_indices(repr_t const& repr);

  dict<IndexBase, type<dolfin::uint> > const index_dimensions_;
  tuple<IndexBase> const indices_;
//      tuple<FixedIndex> const fixed_indices_;

  mutable repr_t repr_;
  mutable std::string str_;

};
} /* namespace ufl */
#endif /* __UFL_INDEX_H_ */
