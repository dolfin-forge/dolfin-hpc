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
    expressions_(fill_expressions(left, right)),
    repr_(*this, left, right),
    str_(left.str() + " " + name + " " + right.str())
  {
  }

//-----------------------------------------------------------------------------
  BinaryCondition::BinaryCondition(std::string const& name, repr_t const & repr):
    Condition(name, repr),
    expressions_(fill_expressions(args())),
    repr_(*this, *expressions_[0], *expressions_[1]),
    str_(expressions_[0]->str() + " " + name + " " + expressions_[1]->str())
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
  std::vector<Expression const *> const BinaryCondition::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr(reprs.size());
    for(dolfin::uint i=0; i<expr.size(); ++i)
      expr[i] = expr[i]->create(reprs[i]);
//    expr.push_back(Expression::create(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const BinaryCondition::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
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
    expressions_(fill_expressions(e)),
    repr_(*this, *expressions_[0]),
    str_("!(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NotCondition::NotCondition(repr_t const & repr):
    Condition("NotCondition", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, *expressions_[0]),
    str_("!(" + expressions_[0]->str() + ")")
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
  std::vector<Expression const *> const NotCondition::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr(reprs.size());
    for(dolfin::uint i=0; i<expr.size(); ++i)
      expr[i] = expr[i]->create(reprs[i]);
//    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NotCondition::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

//-----------------------------------------------------------------------------
  Conditional::Conditional(Condition const& c, Expression const& e1, Expression const& e2) :
    Class("Conditional"),
    c_(&c),
    expressions_(fill_expressions(e1, e2)),
    repr_(*this, *c_, *expressions_[0], *expressions_[1]),
    str_(c_->str() + " ? " + e1.str() + " : " + e2.str())
  {
    std::cout << "Conditional from C++" << std::endl;
    std::cout << repr_ << std::endl;
  }

//-----------------------------------------------------------------------------
  Conditional::Conditional(repr_t const& repr):
    Class("Conditional", repr),
    c_(Condition::create(arg(0))),
    expressions_(fill_expressions(args())),
    repr_(*this, *c_, *expressions_[0], *expressions_[1]),
    str_(c_->str() + " ? " + expressions_[0]->str() + " : " + expressions_[1]->str())
  {
    std::cout << "Conditional from repr" << std::endl;
    std::cout << repr_ << std::endl;
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

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Conditional::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr(2);
    for(dolfin::uint i=0; i<expr.size(); ++i)
      expr[i] = expr[i]->create(reprs[i+1]);
//    expr.push_back(Expression::create(reprs[1]));
//    expr.push_back(Expression::create(reprs[2]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Conditional::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
  }
}
