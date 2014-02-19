// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLConditional.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
  Condition::Condition(std::string const& name) : 
    Class(name)
  {
  }

//-----------------------------------------------------------------------------
  Condition::Condition(std::string const& name, repr_t const & repr) :
    Class(name, repr)
  {
  }

//-----------------------------------------------------------------------------
  Condition::~Condition()
  {
  }

//-----------------------------------------------------------------------------
  Condition const * Condition::create(Object::repr_t const& repr)
  {
    std::string name = Class::make_name(repr);
    if(name == "==")
    {
      return new EQ(repr);
    }
    else if(name == "!=")
    {
      return new NE(repr);
    }
    else if(name == "<=")
    {
      return new LE(repr);
    }
    else if(name == ">=")
    {
      return new GE(repr);
    }
    else if(name == "<")
    {
      return new LT(repr);
    }
    else if(name == ">")
    {
      return new GT(repr);
    }
    else if(name == "&&")
    {
      return new AndCondition(repr);
    }
    else if(name == "||")
    {
      return new OrCondition(repr);
    }
    else if(name == "NotCondition")
    {
      return new NotCondition(repr);
    }
    //is this needed?
//    else if(name == "Conditional")
//    {
//      return new Conditional(repr);
//    }
    else
    {
      error("Unknown type of ufl::Condition: '" + name + "'");
    }
    return NULL;
  }

//-----------------------------------------------------------------------------
//  Object::repr_t const Condition::repr() const
//  {
//    return repr_;
//  }
//
//-----------------------------------------------------------------------------
//  std::string const Condition::str() const
//  {
//    return str_;
//  }
//
//-----------------------------------------------------------------------------
//  void Condition::display() const
//  {
//  }

//-----------------------------------------------------------------------------
  BinaryCondition::BinaryCondition(std::string const& name, Expression const& left, Expression const& right) :
    Condition(name),
    left_expression_(left),
    right_expression_(right),
    repr_(*this, left_expression_, right_expression_),
    str_(left_expression_.str() + " " + name + " " + right_expression_.str())
  {
  }

//-----------------------------------------------------------------------------
  BinaryCondition::BinaryCondition(std::string const& name, repr_t const & repr):
    Condition(name, repr),
    left_expression_(arg(0)),
    right_expression_(arg(1)),
    repr_(*this, left_expression_, right_expression_),
    str_(left_expression_.str() + " " + name + " " + right_expression_.str())
  {
  }

//-----------------------------------------------------------------------------
  BinaryCondition::~BinaryCondition()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const BinaryCondition::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const BinaryCondition::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void BinaryCondition::display() const
  {
  }

//-----------------------------------------------------------------------------
  EQ::EQ(Expression const& e1, Expression const& e2) :
    BinaryCondition("==", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  EQ::EQ(repr_t const& repr) :
    BinaryCondition("==", repr)
  {
  }

//-----------------------------------------------------------------------------
  EQ::~EQ()
  {
  }
  
//-----------------------------------------------------------------------------
  NE::NE(Expression const& e1, Expression const& e2) :
    BinaryCondition("!=", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  NE::NE(repr_t const& repr) :
    BinaryCondition("!=", repr)
  {
  }

//-----------------------------------------------------------------------------
  NE::~NE()
  {
  }
  
//-----------------------------------------------------------------------------
  LE::LE(Expression const& e1, Expression const& e2) :
    BinaryCondition("<=", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  LE::LE(repr_t const & repr):
    BinaryCondition("<=", repr)
  {
  }

//-----------------------------------------------------------------------------
  LE::~LE()
  {
  }
  
//-----------------------------------------------------------------------------
  GE::GE(Expression const& e1, Expression const& e2) :
    BinaryCondition(">=", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  GE::GE(repr_t const & repr):
    BinaryCondition(">=", repr)
  {
  }

//-----------------------------------------------------------------------------
  GE::~GE()
  {
  }
  
//-----------------------------------------------------------------------------
  LT::LT(Expression const& e1, Expression const& e2) :
    BinaryCondition("<", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  LT::LT(repr_t const & repr):
    BinaryCondition("<", repr)
  {
  }

//-----------------------------------------------------------------------------
  LT::~LT()
  {
  }
  
//-----------------------------------------------------------------------------
  GT::GT(Expression const& e1, Expression const& e2) :
    BinaryCondition(">", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  GT::GT(repr_t const & repr):
    BinaryCondition(">", repr)
  {
  }

//-----------------------------------------------------------------------------
  GT::~GT()
  {
  }
  
//-----------------------------------------------------------------------------
  AndCondition::AndCondition(Expression const& e1, Expression const& e2) :
    BinaryCondition("&&", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  AndCondition::AndCondition(repr_t const & repr):
    BinaryCondition("&&", repr)
  {
  }

//-----------------------------------------------------------------------------
  AndCondition::~AndCondition()
  {
  }
  
//-----------------------------------------------------------------------------
  OrCondition::OrCondition(Expression const& e1, Expression const& e2) :
    BinaryCondition("||", e1, e2)
  {
  }

//-----------------------------------------------------------------------------
  OrCondition::OrCondition(repr_t const & repr):
    BinaryCondition("||", repr)
  {
  }

//-----------------------------------------------------------------------------
  OrCondition::~OrCondition()
  {
  }
  
//-----------------------------------------------------------------------------
  NotCondition::NotCondition(Expression const& e) :
    Condition("NotCondition"),
    expression_(e),
    repr_(*this, expression_),
    str_("!(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NotCondition::NotCondition(repr_t const & repr):
    Condition("NotCondition", repr),
    expression_(arg(0)),
    repr_(*this, expression_),
    str_("!(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NotCondition::~NotCondition()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const NotCondition::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const NotCondition::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void NotCondition::display() const
  {
  }

//-----------------------------------------------------------------------------
  Conditional::Conditional(Condition const& c, Expression const& e1, Expression const& e2) :
    Class("Conditional"),
    c_(&c),
    e1_(e1),
    e2_(e2),
    repr_(*this, *c_, e1_, e2_),
    str_(c_->str() + " ? " + e1_.str() + " : " + e2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Conditional::Conditional(repr_t const& repr):
    Class("Conditional", repr),
    c_(Condition::create(arg(0))),
    e1_(arg(1)),
    e2_(arg(2)),
    repr_(*this, *c_, e1_, e2_),
    str_(c_->str() + " ? " + e1_.str() + " : " + e2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Conditional::~Conditional()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Conditional::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Conditional::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Conditional::display() const
  {
  }


}
