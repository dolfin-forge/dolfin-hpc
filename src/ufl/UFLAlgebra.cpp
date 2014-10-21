// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLAlgebra.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Sum::Sum(Expression const& s1, Expression const& s2) :
    Expression("Sum"),
    expressions_(fill_expressions(s1, s2)),
    repr_(*this, expressions_),
    str_(s1.str() + " + " + s2.str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::Sum(repr_t const & repr):
    Expression("Sum", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(expressions_[0]->str() + " + " + expressions_[1]->str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::~Sum()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Sum::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Sum::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0 
      = expressions_[0]->level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands1
      = expressions_[1]->level_operands(operands);

    const dolfin::uint size = std::max(new_operands0.size(), new_operands1.size());
    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      if(i<new_operands0.size())
        for(dolfin::uint j=0; j<new_operands0[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands0[i][j]);
//          std::cout << new_operands0[i][j]->name() << std::endl;
        }
      if(i<new_operands1.size())
        for(dolfin::uint j=0; j<new_operands1[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
        }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  Sum const* Sum::create(Object::repr_t const& repr)
  {
    return new Sum(repr);
  }
  
//-----------------------------------------------------------------------------
  std::vector<Expression const *> const& Sum::operands() const
  {
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Sum::shape() const
  {
    return expressions_[0]->shape();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Sum::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Sum::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Sum::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    expressions_[1]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Sum::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Sum::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Sum::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Sum::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(Expression::create(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Sum::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
  }

//-----------------------------------------------------------------------------
  Product::Product(Expression const& p1, Expression const& p2) :
    Expression("Product"),
    expressions_(fill_expressions(p1, p2)),
    repr_(*this, expressions_),
    str_(p1.str() + " * " + p2.str())
  {
  }

//-----------------------------------------------------------------------------
  Product::Product(repr_t const & repr):
    Expression("Product", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(expressions_[0]->str() + " * " + expressions_[1]->str())
  {
  }

//-----------------------------------------------------------------------------
  Product::~Product()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Product::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Product::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0 
      = expressions_[0]->level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands1
      = expressions_[1]->level_operands(operands);

    const dolfin::uint size = std::max(new_operands0.size(), new_operands1.size());
    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      if(i<new_operands0.size())
        for(dolfin::uint j=0; j<new_operands0[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands0[i][j]);
//          std::cout << new_operands0[i][j]->name() << std::endl;
        }
      if(i<new_operands1.size())
        for(dolfin::uint j=0; j<new_operands1[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
        }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  Product const* Product::create(Object::repr_t const& repr)
  {
    return new Product(repr);
  }
  
//-----------------------------------------------------------------------------
  std::vector<Expression const *> const& Product::operands() const
  {
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Product::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Product::free_indices() const
  {
    std::vector<Index const *> indices;

    for(dolfin::uint i=0; i<expressions_.size(); ++i)
    {
      std::vector<Index const *> const expr_indices = expressions_[i]->free_indices().operands();
      for(dolfin::uint j=0; j<expr_indices.size(); ++j)
        indices.push_back(expr_indices[j]);
    }

    return tuple<Index>(indices);
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Product::index_dimensions() const
  {
    std::vector<std::pair<IndexBase const *, type<dolfin::uint> const *> > i_dims;

    for(dolfin::uint i=0; i<expressions_.size(); ++i)
    {
      std::vector<std::pair<IndexBase const *, type<dolfin::uint> const *> > const expr_idims = expressions_[i]->index_dimensions().map();
      for(dolfin::uint j=0; j<expr_idims.size(); ++j)
        i_dims.push_back(expr_idims[j]);
    }

    return dict<IndexBase, type<dolfin::uint> >(i_dims);
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Product::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::cout << "Product::evaluate" << n << std::endl;
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0 =
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals1 =
    expressions_[1]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > result(new_vals0.size());
    for(dolfin::uint i=0; i<new_vals0.size(); ++i)
    {
      result[i].resize(new_vals0[i].size());
      for(dolfin::uint j=0; j<new_vals0[i].size(); ++j)
      {
        result[i][j].resize(new_vals0[i][j].size());
        for(dolfin::uint k=0; k<new_vals0[i][j].size(); ++k)
          result[i][j][k] = new_vals0[i][j][k] *  new_vals1[i][j][k];
      }
    }

    for(dolfin::uint i=0; i<new_vals0.size(); ++i)
      for(dolfin::uint j=0; j<new_vals0[i].size(); ++j)
        for(dolfin::uint k=0; k<new_vals0[i][j].size(); ++k)
          std::cout << "(" << i << "," << j << "," << k << ") = " << new_vals0[i][j][k] << std::endl; 
    std::cout << std::endl;
    for(dolfin::uint i=0; i<new_vals1.size(); ++i)
      for(dolfin::uint j=0; j<new_vals1[i].size(); ++j)
        for(dolfin::uint k=0; k<new_vals1[i][j].size(); ++k)
          std::cout << "(" << i << "," << j << "," << k << ") = " << new_vals1[i][j][k] << std::endl; 
    std::cout << "Product::evaluate END" << n << std::endl;
    return result;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Product::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Product::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Product::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Product::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(Expression::create(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Product::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
  }

//-----------------------------------------------------------------------------
  Division::Division(Expression const& d1, Expression const& d2) :
    Expression("Division"),
    expressions_(fill_expressions(d1, d2)),
    repr_(*this, expressions_),
    str_(d1.str() + " / " + d2.str())
  {
  }

//-----------------------------------------------------------------------------
  Division::Division(repr_t const & repr):
    Expression("Division", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(expressions_[0]->str() + " / " + expressions_[1]->str())
  {
  }

//-----------------------------------------------------------------------------
  Division::~Division()
  {
  }
  
//-----------------------------------------------------------------------------
  std::vector<Class const*> const Division::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Division::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0 
      = expressions_[0]->level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands1
      = expressions_[1]->level_operands(operands);

    const dolfin::uint size = std::max(new_operands0.size(), new_operands1.size());
    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      if(i<new_operands0.size())
        for(dolfin::uint j=0; j<new_operands0[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands0[i][j]);
//          std::cout << new_operands0[i][j]->name() << std::endl;
        }
      if(i<new_operands1.size())
        for(dolfin::uint j=0; j<new_operands1[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
        }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  Division const* Division::create(Object::repr_t const& repr)
  {
    return new Division(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const& Division::operands() const
  {
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Division::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Division::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Division::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Division::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    expressions_[1]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Division::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Division::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Division::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Division::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(Expression::create(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Division::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
  }

//-----------------------------------------------------------------------------
  Power::Power(Expression const& a, Expression const& b) :
    Expression("Power"),
    expressions_(fill_expressions(a, b)),
    repr_(*this, expressions_),
    str_(a.str() + " ** " + b.str())
  {
  }

//-----------------------------------------------------------------------------
  Power::Power(repr_t const & repr):
    Expression("Power", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(expressions_[0]->str() + " ** " + expressions_[1]->str())
  {
  }

//-----------------------------------------------------------------------------
  Power::~Power()
  {
  }
  
//-----------------------------------------------------------------------------
  std::vector<Class const*> const Power::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Power::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0
      = expressions_[0]->level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands1
      = expressions_[1]->level_operands(operands);

    const dolfin::uint size = std::max(new_operands0.size(), new_operands1.size());
    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      if(i<new_operands0.size())
        for(dolfin::uint j=0; j<new_operands0[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands0[i][j]);
//          std::cout << new_operands0[i][j]->name() << std::endl;
        }
      if(i<new_operands1.size())
        for(dolfin::uint j=0; j<new_operands1[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
        }
    }
    return tmp;
  }

//-----------------------------------------------------------------------------
  Power const* Power::create(Object::repr_t const& repr)
  {
      return new Power(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const& Power::operands() const
  {
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Power::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Power::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Power::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Power::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    expressions_[1]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Power::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Power::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Power::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Power::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(Expression::create(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Power::fill_expressions(Expression const& e1, Expression const& e2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    expr.push_back(&e2);
    return expr;
  }

//-----------------------------------------------------------------------------
  Abs::Abs(Expression const& a) :
    Expression("Abs"),
    expressions_(fill_expressions(a)),
    repr_(*this, expressions_),
    str_("| " + a.str() + " |")
  {
  }

//-----------------------------------------------------------------------------
  Abs::Abs(repr_t const & repr):
    Expression("Abs", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("| " + expressions_[0]->str() + " |")
  {
  }

//-----------------------------------------------------------------------------
  Abs::~Abs()
  {
  }
  
//-----------------------------------------------------------------------------
  std::vector<Class const*> const Abs::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Abs::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0 
      = expressions_[0]->level_operands(operands);

    const dolfin::uint size = new_operands0.size();
    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      for(dolfin::uint j=0; j<new_operands0[i].size(); ++j)
      {
        tmp[i+1].push_back(new_operands0[i][j]);
//        std::cout << new_operands0[i][j]->name() << std::endl;
      }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  Abs const* Abs::create(Object::repr_t const& repr)
  {
      return new Abs(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const& Abs::operands() const
  {
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Abs::shape() const
  {
    return expressions_[0]->shape();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Abs::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Abs::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Abs::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Abs::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Abs::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Abs::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Abs::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Abs::fill_expressions(Expression const& e1)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e1);
    return expr;
  }

}
