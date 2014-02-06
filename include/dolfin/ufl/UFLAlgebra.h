// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_ALGEBRA_H_
#define __UFL_ALGEBRA_H_

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
//#include <dolfin/ufl/UFLIndex.h>

//#include <dolfin/common/types.h>

namespace ufl
{

  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Product.
   */

  class Product : public Class
  {

    public:

      ///
      Product(Expression const& p1, Expression const& p2);

//      ///
//      Product(std::vector<Expression const> const& p);

      ///
      ~Product();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Expression const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      Product const* create(repr_t const & repr) const;

    private:

      Expression const p1_;
      Expression const p2_;

      mutable repr_t repr_;
      mutable std::string str_;

  };
} /* namespace ufl */
#endif /* __UFL_ALGEBRA_H_ */
