// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLAlgebra.h>
#include <dolfin/ufl/UFLDifferentiation.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndexed.h>
#include <dolfin/ufl/UFLIndexSum.h>
#include <dolfin/ufl/UFLTensors.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Expression::Expression(Object const& object) :
    is_cellwise_constant_(false)
  {
//    std::stringstream ssrepr;
//    repr_ = ssrepr.str();
//
//    std::stringstream ssstr;
//    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  Expression::Expression(repr_t const & repr) :
    is_cellwise_constant_(false)
  {
//    if(repr.length() == 0)
//      dolfin_assert("An empty signature was passed to create an Expression.");
//
//
//    Object const * object;
//
//    std::string::const_iterator it;
//    std::string help_string;
//    dolfin::uint i = 0;
//    for(it = repr.begin(); it!=repr.end(); ++it, ++i)
//    {
//      help_string += *it;
//
//      if(help_string == "ComponentTensor")
//      {
//        std::cout << "create ComponentTensor form Expression : " << repr << std::endl;
//        ComponentTensor const * c_tensor;
//        c_tensor = c_tensor->create(repr);
//        object = c_tensor;
//        return new Expression(*object);
//      }
//
//      if(help_string == "IndexSum")
//      {
//        std::cout << "create IndexSum form Expression : " << repr << std::endl;
//        IndexSum const * index_sum;
//        index_sum = index_sum->create(repr);
//        object = index_sum;
//        return new Expression(*object);
//      }
//
//      if(help_string == "Indexed")
//      {
//        std::cout << "create Indexed form Expression : " << repr << std::endl;
//        Indexed const * indexed;
//        indexed = indexed->create(repr);
//        object = indexed;
//        return new Expression(*object);
//      }
//
//      if(help_string == "Product")
//      {
//        std::cout << "create Product form Expression : " << repr << std::endl;
//        Product const * p;
//        p = p->create(repr);
//        object = p;
//        return new Expression(*object);
//      }
//
//      if(help_string == "SpatialDerivative")
//      {
//        std::cout << "create SpatialDerivative form Expression : " << repr << std::endl;
//        SpatialDerivative const * sd;
//        sd = sd->create(repr);
//        object = sd;
//        return new Expression(*object);
//      }
//    }
//
//    return new Expression(*object);
  }
//-----------------------------------------------------------------------------
  Expression::~Expression()
  {
  }
  
//-----------------------------------------------------------------------------
  Cell const& Expression::cell() const
  {
  }

//-----------------------------------------------------------------------------
  bool const Expression::is_cellwise_constant() const
  {
    return is_cellwise_constant_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Expression::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const Expression::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Expression::display() const
  {
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
