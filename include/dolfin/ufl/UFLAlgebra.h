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
   *  @brief  Provides an interface complying with UFL Sum.
   */

  class Sum : public Class
  {

    public:

      ///
      Sum(Expression const& s1, Expression const& s2);

      ///
      Sum(repr_t const& repr);

      ///
      ~Sum();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Expression const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const s1_;
      Expression const s2_;

      repr_t const repr_;
      std::string const str_;
  };


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

      ///
      Product (repr_t const & repr);

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

      Expression const p1_;
      Expression const p2_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Division.
   */

  class Division : public Class
  {

    public:

      ///
      Division(Expression const& d1, Expression const& d2);

      ///
      Division (repr_t const & repr);

      ///
      ~Division();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Expression const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const d1_;
      Expression const d2_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Power.
   */

  class Power : public Class
  {

    public:

      ///
      Power(Expression const& a, Expression const& b);

      ///
      Power (repr_t const & repr);

      ///
      ~Power();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Expression const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const a_;
      Expression const b_;

      repr_t const repr_;
      std::string const str_;
  };

  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Abs.
   */

  class Abs : public Class
  {

    public:

      ///
      Abs(Expression const& a);

      ///
      Abs (repr_t const & repr);

      ///
      ~Abs();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const a_;

      repr_t const repr_;
      std::string const str_;
  };
} /* namespace ufl */
#endif /* __UFL_ALGEBRA_H_ */
