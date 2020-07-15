// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLAlgebra.h>
#include <dolfin/ufl/UFLArgument.h>
#include <dolfin/ufl/UFLCoefficient.h>
#include <dolfin/ufl/UFLDifferentiation.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndexed.h>
#include <dolfin/ufl/UFLIndexSum.h>
#include <dolfin/ufl/UFLTensors.h>

#include <dolfin/log/log.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
Expression::Expression(std::string const& name) :
    Class(name),
    is_cellwise_constant_(true)
{
}

//-----------------------------------------------------------------------------
Expression::Expression(std::string const& name, repr_t const & repr) :
    Class(name, repr),
    is_cellwise_constant_(true)
{
}

//-----------------------------------------------------------------------------
Expression::~Expression()
{
}

//-----------------------------------------------------------------------------
Cell const Expression::cell() const
{
  for (dolfin::uint i = 0; i < operands().size(); ++i)
  {
    Cell const& d = operands()[i]->cell();
    if (d.domain().type() != Domain::None) return d;
  }
  
  return Cell(Domain::None);
}

//-----------------------------------------------------------------------------
dolfin::uint Expression::geometric_dimension() const
{
  return cell().geometric_dimension();
}

//-----------------------------------------------------------------------------
bool Expression::is_cellwise_constant() const
{
  return is_cellwise_constant_;
}

//-----------------------------------------------------------------------------
Expression const* Expression::create(Object::repr_t const& repr)
{
  std::string name = Class::make_name(repr);
  if (name == "Sum")
  {
    return new Sum(repr);
  }
  else if (name == "Product")
  {
    return new Product(repr);
  }
  else if (name == "Division")
  {
    return new Division(repr);
  }
  else if (name == "Power")
  {
    return new Power(repr);
  }
  else if (name == "Abs")
  {
    return new Abs(repr);
  }
  else if (name == "Argument")
  {
    return new Argument(repr);
  }
  else if (name == "Coefficient")
  {
    return new Coefficient(repr);
  }
  else if (name == "Constant")
  {
    return new Constant(repr);
  }
  else if (name == "VectorConstant")
  {
    return new VectorConstant(repr);
  }
  else if (name == "TensorConstant")
  {
    return new TensorConstant(repr);
  }
  else if (name == "CoefficientDerivative")
  {
    return new CoefficientDerivative(repr);
  }
  else if (name == "SpatialDerivative")
  {
    return new SpatialDerivative(repr);
  }
  else if (name == "VariableDerivative")
  {
    return new VariableDerivative(repr);
  }
  else if (name == "Grad")
  {
    return new Grad(repr);
  }
  else if (name == "Div")
  {
    return new Div(repr);
  }
  else if (name == "NablaGrad")
  {
    return new NablaGrad(repr);
  }
  else if (name == "NablaDiv")
  {
    return new NablaDiv(repr);
  }
  else if (name == "Curl")
  {
    return new Curl(repr);
  }
  else if (name == "Index")
  {
    return new Index(repr);
  }
  else if (name == "FixedIndex")
  {
    return new FixedIndex(repr);
  }
  else if (name == "MultiIndex")
  {
    return new MultiIndex(repr);
  }
  else if (name == "IndexSum")
  {
    return new IndexSum(repr);
  }
  else if (name == "Indexed")
  {
    return new Indexed(repr);
  }
  else if (name == "ComponentTensor")
  {
    return new ComponentTensor(repr);
  }
  else if (name == "Tuple")
  {
    return new Tuple(repr);
  }
  else if (name == "Label")
  {
    return new Label(repr);
  }
  else if (name == "Variable")
  {
    return new Variable(repr);
  }
  else if (name == "Data")
  {
    return new Data(repr);
  }
  else
  {
    error("Unknown type of ufl::Expression: '" + name + "'");
  }
  
  return nullptr;
}

/*
 //-----------------------------------------------------------------------------
 Operator::Operator(Expression const& expression) :
 expression_(expression)
 {
 }

 //-----------------------------------------------------------------------------
 Operator::~Operator()
 {
 }

 //-----------------------------------------------------------------------------
 bool const Operator::is_cellwise_constant() const
 {
 return expression_.is_cellwise_constant();
 }
 */
}
