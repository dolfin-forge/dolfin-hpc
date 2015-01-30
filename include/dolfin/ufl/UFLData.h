// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_DATA_H_
#define __UFL_DATA_H_

#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Data
 *
 *  @brief  Provides an interface complying with UFL Data.
 */

  class Data : public Expression
  {
    public:

      ///
      Data(Expression const& expression);

      ///
      Data(repr_t const& repr);

      ///
      ~Data();

      ///
      virtual std::vector<Class const* > const operands (std::string const& name) const;

      ///
      virtual std::vector<std::vector<Class const *> > const level_operands (
          std::vector<std::vector<Class const *> > const& operands) const;

      //--- INTERFACE -------------------------------------------------------------

      ///
      static Data const * create (Object::repr_t const& repr);

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
          ufc::cell const& ref_cell, 
          std::vector<dolfin::real*> const& q_points,
          const double * const * coordinates) const; 

      //--- INTERFACE inherited from UFLClass -------------------------------------
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

    private:

      std::vector<Expression const *> const fill_expressions(std::vector<repr_t> const& reprs);
      std::vector<Expression const *> const expressions_;

      repr_t const repr_;
      std::string const str_;

  };

} /* namespace ufl */
#endif /* __UFL_DATA_H_ */
