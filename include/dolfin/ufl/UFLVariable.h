// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_VARIABLE_H_
#define __UFL_VARIABLE_H_

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLtype.h>

namespace ufl
{

  /**
   *  DOCUMENTATION:
   *
   *  @class  Label
   *
   *  @brief  Provides an interface complying with UFL Label.
   */

  class Label : public Class
  {
    public:

      ///
      Label(dolfin::uint const& count);

      ///
      Label (repr_t const & repr);

      ///
      ~Label();

      //--- INTERFACE -------------------------------------------------------------

      type<dolfin::uint> const& count() const;
      
      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      type<dolfin::uint> count_;

      repr_t const repr_;
      std::string const str_;
  };

  /**
   *  DOCUMENTATION:
   *
   *  @class  Variable
   *
   *  @brief  Provides an interface complying with UFL Variable.
   */

  class Variable : public Class
  {
    public:

      ///
      Variable(Expression const& expression, Label const& label);

      ///
      Variable (repr_t const & repr);

      ///
      ~Variable();

      //--- INTERFACE -------------------------------------------------------------

      std::pair<Expression, Label> const& operands() const;
      
//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

      /// UFL: Return whether this expression is spatially constant over each cell
      bool const is_cellwise_constant() const;

//      ///
//      evaluate(self, x, mapping, component, index_values):
      
      
      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      std::pair<Expression, Label> const expr_label_;

      repr_t const repr_;
      std::string const str_;
  };
} /* namespace ufl */
#endif /* __UFL_VARIABLE_H_ */
