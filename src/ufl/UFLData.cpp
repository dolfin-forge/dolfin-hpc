// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLData.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
Data::Data(Expression const& expr) :
    Expression("Data"),
    expressions_(std::vector<Expression const *>(1, &expr)),
    repr_(*this, expressions_),
    str_("Data(" + expr.str() + ")")
{
}

//-----------------------------------------------------------------------------
Data::Data(repr_t const& repr) :
    Expression("Data", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("Data(" + expressions_[0]->str() + ")")
{
}

//-----------------------------------------------------------------------------
Data::~Data()
= default;

//-----------------------------------------------------------------------------
std::vector<Class const*> const Data::operands(std::string const& name) const
{
  std::vector<Class const*> expr0 = expressions_[0]->operands(name);
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Data::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0 =
      expressions_[0]->level_operands(operands);
  
  const dolfin::uint size = new_operands0.size();
  std::vector<std::vector<Class const*> > tmp(size + 1);
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  
  tmp[0] = obj0;
  for (dolfin::uint i = 0; i < tmp.size() - 1; ++i)
  {
    for (dolfin::uint j = 0; j < new_operands0[i].size(); ++j)
    {
      tmp[i + 1].push_back(new_operands0[i][j]);
//        std::cout << new_operands0[i][j]->name() << std::endl;
    }
  }
  
  return tmp;
}

//-----------------------------------------------------------------------------
Data const* Data::create(Object::repr_t const& repr)
{
  return new Data(repr);
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Data::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
ValueArray const Data::shape() const
{
  return ValueArray();
}

//-----------------------------------------------------------------------------
tuple<Index> const Data::free_indices() const
{
  return tuple<Index>();
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const Data::index_dimensions() const
{
  return dict<IndexBase, type<dolfin::uint> >();
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const Data::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Data::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Data::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Data::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Data::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  expr.push_back(Expression::create(reprs[0]));
  return expr;
}
}
