// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLDifferentiation.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  CoefficientDerivative::CoefficientDerivative(type<dolfin::real> const& integrand,
      Tuple const& coefficients, Tuple const& arguments/*, 
      type<std::map<dolfin::uint, dolfin::uint> > const& coeff_derivatives*/) :
    Class("CoefficientDerivative"),
    integrand_(integrand),
    coefficients_(coefficients),
    arguments_(arguments),
//    coeff_derivatives_(coeff_derivatives),
    repr_(*this, integrand_, coefficients_, arguments_/*, coeff_derivatives_*/),
    str_("d/dfj { " + integrand_.str() + " } with fh=" + coefficients_.str() + ", dfh/dfj = " + arguments_.str() + ". and coefficient derivatives "/* + coeff_derivatives_.str()*/)
  {
  }

//-----------------------------------------------------------------------------
  CoefficientDerivative::CoefficientDerivative(repr_t const & repr) :
    Class("CoefficientDerivative"),
    integrand_(arg(0)),
    coefficients_(arg(1)),
    arguments_(arg(2)),
//    coeff_derivatives_(arg(3)),
    repr_(*this, integrand_, coefficients_, arguments_/*, coeff_derivatives_*/),
    str_("d/dfj { " + integrand_.str() + " } with fh=" + coefficients_.str() + ", dfh/dfj = "
        + arguments_.str() + ". and coefficient derivatives "/* + coeff_derivatives_.str()*/)
  {
  }

//-----------------------------------------------------------------------------
  CoefficientDerivative::~CoefficientDerivative()
  {
  }
  
//-----------------------------------------------------------------------------
//  std::pair<Expression const, Index const> const& CoefficientDerivative::operands() const
//  {
//    return std::make_pair(expression_, index_);  
//  }

//-----------------------------------------------------------------------------
  Object::repr_t const CoefficientDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const CoefficientDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void CoefficientDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::SpatialDerivative(Expression const& expression, Index const& index) :
    Class("SpatialDerivative"),
    expr_index_(expression, index),
    repr_(*this, expr_index_.first, expr_index_.second),
    str_("d/dx_" + expr_index_.second.str() + " " + expr_index_.first.str())
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::SpatialDerivative(repr_t const & repr) :
    Class("SpatialDerivative"),
    expr_index_(arg(0), arg(1)),
    repr_(*this, expr_index_.first, expr_index_.second),
    str_("d/dx_" + expr_index_.second.str() + " " + expr_index_.first.str())
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::~SpatialDerivative()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Index> const& SpatialDerivative::operands() const
  {
    return expr_index_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const SpatialDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const SpatialDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void SpatialDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  VariableDerivative::VariableDerivative(Expression const& expression, Variable const& variable) :
    Class("VariableDerivative"),
    expr_var_(expression, variable),
    repr_(*this, expr_var_.first, expr_var_.second),
    str_("d/d[" + expr_var_.second.str() + "] " + expr_var_.first.str())
  {
  }

//-----------------------------------------------------------------------------
  VariableDerivative::VariableDerivative(repr_t const & repr) :
    Class("VariableDerivative"),
    expr_var_(arg(0), arg(1)),
    repr_(*this, expr_var_.first, expr_var_.second),
    str_("d/d[" + expr_var_.second.str() + "] " + expr_var_.first.str())
  {
  }

//-----------------------------------------------------------------------------
  VariableDerivative::~VariableDerivative()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Variable> const& VariableDerivative::operands() const
  {
    return expr_var_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const VariableDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const VariableDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void VariableDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  Grad::Grad(Expression const& expression) :
    Class("Grad"),
    expression_(expression),
    repr_(*this, expression_),
    str_("grad(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Grad::Grad(repr_t const & repr) :
    Class("Grad"),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("grad(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Grad::~Grad()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& Grad::operands() const
  {
    return expression_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Grad::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Grad::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Grad::display() const
  {
  }

//-----------------------------------------------------------------------------
  Div::Div(Expression const& expression) :
    Class("Div"),
    expression_(expression),
    repr_(*this, expression_),
    str_("div(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Div::Div(repr_t const & repr) :
    Class("Div"),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("div(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Div::~Div()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& Div::operands() const
  {
    return expression_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Div::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Div::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Div::display() const
  {
  }

//-----------------------------------------------------------------------------
  NablaGrad::NablaGrad(Expression const& expression) :
    Class("NablaGrad"),
    expression_(expression),
    repr_(*this, expression_),
    str_("nabla_grad(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaGrad::NablaGrad(repr_t const & repr) :
    Class("NablaGrad"),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("nabla_grad(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaGrad::~NablaGrad()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& NablaGrad::operands() const
  {
    return expression_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const NablaGrad::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const NablaGrad::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void NablaGrad::display() const
  {
  }

//-----------------------------------------------------------------------------
  NablaDiv::NablaDiv(Expression const& expression) :
    Class("NablaDiv"),
    expression_(expression),
    repr_(*this, expression_),
    str_("nabla_div(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaDiv::NablaDiv(repr_t const & repr) :
    Class("NablaDiv"),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("nabla_div(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaDiv::~NablaDiv()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& NablaDiv::operands() const
  {
    return expression_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const NablaDiv::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const NablaDiv::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void NablaDiv::display() const
  {
  }

//-----------------------------------------------------------------------------
  Curl::Curl(Expression const& expression) :
    Class("Curl"),
    expression_(expression),
    repr_(*this, expression_),
    str_("curl(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Curl::Curl(repr_t const & repr) :
    Class("Curl"),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("curl(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Curl::~Curl()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& Curl::operands() const
  {
    return expression_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Curl::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Curl::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Curl::display() const
  {
  }
}
