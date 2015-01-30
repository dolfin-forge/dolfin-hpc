// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIndexed.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Indexed::Indexed(Expression const& expression, MultiIndex const& multi_index) :
    Expression("Indexed"),
    expressions_(fill_expressions(expression, multi_index)),
    repr_(*this, expressions_),
    str_(expression.str() + "[" + multi_index.str() + "]")
  {
  }

//-----------------------------------------------------------------------------
  Indexed::Indexed(Expression const& expression, tuple<IndexBase> const& indices,
      dict<IndexBase, type<dolfin::uint> > const& index_dimensions) :
    Expression("Indexed"),
    expressions_(fill_expressions(expression, indices, index_dimensions)),
    repr_(*this, expressions_),
    str_(expression.str() + "[" + expressions_[1]->str() + "]")
  {
  }

//-----------------------------------------------------------------------------
//  Indexed::Indexed(Expression const& expression, tuple<FixedIndex> const& indices) :
//    Expression("Indexed"),
//    expressions_(fill_expressions(expression, indices)),
//    repr_(*this, expressions_),
//    str_(expression.str() + "[" + expressions_[1]->str() + "]")
//  {
//  }

//-----------------------------------------------------------------------------
  Indexed::Indexed(repr_t const& repr) :
    Expression("Indexed", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(expressions_[0]->str() + "[" + expressions_[1]->str() + "]")
  {
  }

//-----------------------------------------------------------------------------
  Indexed::~Indexed()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Indexed::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Indexed::level_operands(
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
  Indexed const* Indexed::create(Object::repr_t const& repr)
  {
    return new Indexed(repr);
  }
    
//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Indexed::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const Indexed::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Indexed::free_indices() const
  {
    std::vector<Index const *> new_indices = 
      expressions_[0]->free_indices().operands();

    std::vector<Index const *> const& f_indices = 
      expressions_[1]->free_indices().operands();

    for(dolfin::uint i=0; i<f_indices.size(); ++i)
        new_indices.push_back(f_indices[i]);

    return tuple<Index>(new_indices);
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Indexed::index_dimensions() const
  {
    ValueArray expr_shape = expressions_[0]->shape();
    std::vector<Index const *> const& free_indices =
      expressions_[1]->free_indices().operands();

    if(expr_shape.size() != free_indices.size())
      error("Invalid number of indices (%d) for t"\
          "ensor expression of rank %d", free_indices.size(), expr_shape.size());
    
    std::vector<std::pair<IndexBase const *, type<dolfin::uint> const*> > map;

    for(dolfin::uint i=0; i<free_indices.size(); ++i)
      for(dolfin::uint s=0; s<expr_shape.size(); ++s)
	map.push_back(std::pair<IndexBase const *, type<dolfin::uint> const *>(free_indices[i], new type<dolfin::uint>(expr_shape[s])));

    return dict<IndexBase, type<dolfin::uint> >(map);
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Indexed::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::cout << "Indexed::evaluate" << n << std::endl;
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0 =
    expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    for(dolfin::uint i=0; i<new_vals0.size(); ++i)
      for(dolfin::uint j=0; j<new_vals0[i].size(); ++j)
        for(dolfin::uint k=0; k<new_vals0[i][j].size(); ++k)
          std::cout << "(" << i << "," << j << "," << k << ") = " << new_vals0[i][j][k] << std::endl; 
    std::cout << "Indexed::evaluate END" << n << std::endl;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Indexed::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Indexed::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Indexed::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Indexed::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(new MultiIndex(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Indexed::fill_expressions(Expression const& expression, MultiIndex const& multi_index)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&expression);
    expr.push_back(&multi_index);
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Indexed::fill_expressions(Expression const& expression,
      tuple<IndexBase> const& indices, dict<IndexBase, type<dolfin::uint> > const& index_dimensions)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&expression);
    expr.push_back(new MultiIndex(indices, index_dimensions));
    return expr;
  }

//-----------------------------------------------------------------------------
//  std::vector<Expression const *> const Indexed::fill_expressions(Expression const& expression,
//      tuple<FixedIndex> const& indices)
//  {
//    std::vector<Expression const *> expr;
//    expr.push_back(&expression);
//    expr.push_back(new MultiIndex(indices));
//    return expr;
//  }

}
