// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_DIFFERENTIATION_H_
#define __UFL_DIFFERENTIATION_H_

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>
#include <dolfin/ufl/UFLTuple.h>
#include <dolfin/ufl/UFLVariable.h>

//#include <dolfin/common/types.h>

#include <map>

namespace ufl
{

  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL CoefficientDerivative.
   */

  class CoefficientDerivative : public Class
  {

    public:

      ///
      CoefficientDerivative(type<dolfin::real> const& integrand, Tuple const& coefficients, 
          Tuple const& arguments/*, type<std::map<dolfin::uint, dolfin::uint> > const& coeff_derivatives*/);

      ///
      CoefficientDerivative(repr_t const& repr);

      ///
      ~CoefficientDerivative();

      //--- INTERFACE -------------------------------------------------------------

      ///
//      std::pair<Expression const, Index const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      type<dolfin::real> const integrand_;
      Tuple const coefficients_; 
      Tuple const arguments_;
//      std::map<dolfin::uint, dolfin::uint> const& coeff_derivatives_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL SpatialDerivative.
   */

  class SpatialDerivative : public Class
  {

    public:

      ///
      SpatialDerivative(Expression const& expression, Index const& index);

      ///
      SpatialDerivative(repr_t const& repr);

      ///
      ~SpatialDerivative();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Index const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

//      ///
//      evaluate(self, x, mapping, component, index_values, derivatives=()):

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

      Expression const expression_;
      Index const index_;

      repr_t const repr_;
      std::string const str_;

  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL VariableDerivative.
   */

  class VariableDerivative : public Class
  {

    public:

      ///
      VariableDerivative(Expression const& expression, Variable const& variable);

      ///
      VariableDerivative(repr_t const& repr);

      ///
      ~VariableDerivative();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, Variable const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;
      Variable const variable_;

      repr_t const repr_;
      std::string const str_;

  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Grad.
   */

  class Grad : public Class
  {

    public:

      ///
      Grad(Expression const& expression);

      ///
      Grad(repr_t const& repr);

      ///
      ~Grad();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Div.
   */

  class Div : public Class
  {

    public:

      ///
      Div(Expression const& expression);

      ///
      Div(repr_t const& repr);

      ///
      ~Div();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;

      repr_t const repr_;
      std::string const str_;

  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL NablaGrad.
   */

  class NablaGrad : public Class
  {

    public:

      ///
      NablaGrad(Expression const& expression);

      ///
      NablaGrad(repr_t const& repr);

      ///
      ~NablaGrad();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL NablaDiv.
   */

  class NablaDiv : public Class
  {

    public:

      ///
      NablaDiv(Expression const& expression);

      ///
      NablaDiv(repr_t const& repr);

      ///
      ~NablaDiv();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;

      repr_t const repr_;
      std::string const str_;
  };


  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Curl.
   */

  class Curl : public Class
  {

    public:

      ///
      Curl(Expression const& expression);

      ///
      Curl(repr_t const& repr);

      ///
      ~Curl();

      //--- INTERFACE -------------------------------------------------------------

      ///
      Expression const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

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

      Expression const expression_;

      repr_t const repr_;
      std::string const str_;
  };
} /* namespace ufl */
#endif /* __UFL_DIFFERENTIATION_H_ */
