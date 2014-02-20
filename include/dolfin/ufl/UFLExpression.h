// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_EXPRESSION_H_
#define __UFL_EXPRESSION_H_

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLIntegral.h>

namespace ufl
{

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
      Expression(Object const& object);

      ///
      Expression (repr_t const & repr);

//      ///
//      Expression();
      
      ///
      ~Expression();

      //--- INTERFACE -------------------------------------------------------------

//      std::vector<Expression> operands() const;
      
      /// UFL: Return the tensor shape of the expression
//      shape();
      
      /// UFL: Return the tensor rank of the expression
//      uint const rank() const;
      
      /// UFL: Return the cell this expression is defined on
      Cell const& cell() const;
      
      /// UFL: Return the geometric dimension this expression is defined on
//      uint const geometric_dimension() const;
      
      /// UFL: Return whether this expression is spatially constant over each cell
      bool const is_cellwise_constant() const;
      
      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      repr_t const repr_;
      std::string const str_;

      bool const is_cellwise_constant_;
  };

  /*
  class Operator : public Expression
  {
    public:

      Operator (Expression const& expression);
      ~Operator ();

      bool const is_cellwise_constant() const;
    
    private:
      Expression const expression_;
  };
  */
} /* namespace ufl */
#endif /* __UFL_EXPRESSION_H_ */
