// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLTensors.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
ComponentTensor::ComponentTensor(Expression const& expression,
                                 MultiIndex const& multi_index) :
    Expression("ComponentTensor"),
    expressions_(fill_expressions(expression, multi_index)),
    repr_(*this, expressions_),
    str_("{ A | A_{" + multi_index.str() + "} = " + expression.str() + " }")
{
}

//-----------------------------------------------------------------------------
ComponentTensor::ComponentTensor(
    Expression const& expression, tuple<IndexBase> const& indices,
    dict<IndexBase, type<dolfin::uint> > const& index_dimensions) :
    Expression("ComponentTensor"),
    expressions_(fill_expressions(expression, indices, index_dimensions)),
    repr_(*this, expressions_),
    str_(
        "{ A | A_{" + expressions_[1]->str() + "} = " + expressions_[0]->str()
            + " }")
{
}

//-----------------------------------------------------------------------------
//  ComponentTensor::ComponentTensor(Expression const& expression, tuple<FixedIndex> const& indices) :
//    Expression("ComponentTensor"),
//    expressions_(fill_expressions(expression, indices)),
//    repr_(*this, expressions_),
//    str_("{ A | A_{" + expressions_[1]->str() + "} = " + expressions_[0]->str() + " }")
//  {
//  }

//-----------------------------------------------------------------------------
ComponentTensor::ComponentTensor(repr_t const & repr) :
    Expression("ComponentTensor", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_(
        "{ A | A_{" + expressions_[1]->str() + "} = " + expressions_[0]->str()
            + " }")
{
}

//-----------------------------------------------------------------------------
ComponentTensor::~ComponentTensor()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const ComponentTensor::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0 = expressions_[0]->operands(name);
  std::vector<Class const*> expr1 = expressions_[1]->operands(name);
  expr0.insert(expr0.end(), expr1.begin(), expr1.end());
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const ComponentTensor::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0 =
      expressions_[0]->level_operands(operands);
  std::vector<std::vector<Class const*> > new_operands1 =
      expressions_[1]->level_operands(operands);
  
  const dolfin::uint size = std::max(new_operands0.size(),
                                     new_operands1.size());
  std::vector<std::vector<Class const*> > tmp(size + 1);
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  
  tmp[0] = obj0;
  for (dolfin::uint i = 0; i < tmp.size() - 1; ++i)
  {
    if (i < new_operands0.size()) for (dolfin::uint j = 0;
        j < new_operands0[i].size(); ++j)
    {
      tmp[i + 1].push_back(new_operands0[i][j]);
//          std::cout << new_operands0[i][j]->name() << std::endl;
    }
    if (i < new_operands1.size()) for (dolfin::uint j = 0;
        j < new_operands1[i].size(); ++j)
    {
      tmp[i + 1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
    }
  }
  
  return tmp;
}

//-----------------------------------------------------------------------------
ComponentTensor const* ComponentTensor::create(Object::repr_t const& repr)
{
  return new ComponentTensor(repr);
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const ComponentTensor::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
ValueArray const ComponentTensor::shape() const
{
  ValueArray shape_array;
  
  std::vector<std::pair<IndexBase const *, type<dolfin::uint> const *> >
  const& dict_map = expressions_[1]->index_dimensions().map();
  for (dolfin::uint i = 0; i < dict_map.size(); ++i)
    shape_array.push_back(*dict_map[i].second);
  return shape_array;
}

//-----------------------------------------------------------------------------
tuple<Index> const ComponentTensor::free_indices() const
{
  std::vector<Index const *> const& expr_indices =
      expressions_[0]->free_indices().operands();
  
  std::vector<Index const *> const& self_indices =
      expressions_[1]->free_indices().operands();
  
  std::vector<Index const *> new_indices;
  for (dolfin::uint i = 0; i < expr_indices.size(); ++i)
    for (dolfin::uint j = 0; j < self_indices.size(); ++j)
      if (expr_indices[i] != self_indices[j]) new_indices.push_back(
          expr_indices[i]);
  
  return tuple<Index>(new_indices);
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const ComponentTensor::index_dimensions() const
{
  std::vector<Index const *> const& f_indices = free_indices().operands();
  std::vector<std::pair<IndexBase const *, type<dolfin::uint> const*> > const& dims =
      expressions_[0]->index_dimensions().map();
  
  if (f_indices.size() != dims.size()) error(
      "Size of free indices does not match size of index dimensions.");
  
  std::vector<std::pair<IndexBase const *, type<dolfin::uint> const*> > map;
  
  for (dolfin::uint i = 0; i < f_indices.size(); ++i)
    map.push_back(
        std::pair<IndexBase const *, type<dolfin::uint> const*>(
            f_indices[i], dims[i].second));
  
  return dict<IndexBase, type<dolfin::uint> >(map);
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const ComponentTensor::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  std::cout << "ComponentTensor::evaluate" << n << std::endl;
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0 =
      expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
  for (dolfin::uint i = 0; i < new_vals0.size(); ++i)
    for (dolfin::uint j = 0; j < new_vals0[i].size(); ++j)
      for (dolfin::uint k = 0; k < new_vals0[i][j].size(); ++k)
        std::cout << "(" << i << "," << j << "," << k << ") = "
                  << new_vals0[i][j][k] << std::endl;
  std::cout << "ComponentTensor::evaluate END" << n << std::endl;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Object::repr_t const& ComponentTensor::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& ComponentTensor::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void ComponentTensor::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const ComponentTensor::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  expr.push_back(Expression::create(reprs[0]));
  expr.push_back(new MultiIndex(reprs[1]));
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const ComponentTensor::fill_expressions(
    Expression const& e, MultiIndex const& i)
{
  std::vector<Expression const *> expr;
  expr.push_back(&e);
  expr.push_back(&i);
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const ComponentTensor::fill_expressions(
    Expression const& expression, tuple<IndexBase> const& indices,
    dict<IndexBase, type<dolfin::uint> > const& index_dimensions)
{
  std::vector<Expression const *> expr;
  expr.push_back(&expression);
  expr.push_back(new MultiIndex(indices, index_dimensions));
  return expr;
}

//-----------------------------------------------------------------------------
//  std::vector<Expression const *> const ComponentTensor::fill_expressions(Expression const& expression,
//      tuple<FixedIndex> const& indices)
//  {
//    std::vector<Expression const *> expr;
//    expr.push_back(&expression);
//    expr.push_back(new MultiIndex(indices));
//    return expr;
//  }

}
