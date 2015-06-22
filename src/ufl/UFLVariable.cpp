// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/ufl/UFLVariable.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Label::Label(dolfin::uint const& count) :
    Expression("Label"),
    count_(count),
    repr_(*this, count_),
    str_("Label(" + count_.str() + ")")
{
}

//-----------------------------------------------------------------------------
Label::Label(repr_t const& repr) :
    Expression("Label", repr),
    count_(arg(0)),
    repr_(*this, count_),
    str_("Label(" + count_.str() + ")")
{
}
//-----------------------------------------------------------------------------
Label::~Label()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Label::operands(std::string const& name) const
{
  std::vector<Class const*> expr0;
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Label::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);
  
  for (dolfin::uint i = 0; i < operands.size(); ++i)
    new_operands0.push_back(operands[i]);
  
  return new_operands0;
}

//-----------------------------------------------------------------------------
Label const* Label::create(Object::repr_t const& repr) const
{
  return new Label(repr);
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Label::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
ValueArray const Label::shape() const
{
  return ValueArray();
}

//-----------------------------------------------------------------------------
tuple<Index> const Label::free_indices() const
{
  return tuple<Index>();
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const Label::index_dimensions() const
{
  return dict<IndexBase, type<dolfin::uint> >();
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const Label::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& Label::count() const
{
  return count_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Label::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Label::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Label::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Label::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Label::fill_expressions(
    dolfin::uint const& count)
{
  std::vector<Expression const *> expr;
  return expr;
}

//-----------------------------------------------------------------------------
Variable::Variable(Expression const& expression, Label const& label) :
    Expression("Variable"),
    expressions_(fill_expressions(expression, label)),
//    label_(label),
    repr_(*this, expression, label),
    str_("var" + label.count().str() + "(" + expressions_[0]->str() + ")")
{
}

//-----------------------------------------------------------------------------
Variable::Variable(repr_t const& repr) :
    Expression("Variable", repr),
    expressions_(fill_expressions(args())),
//    label_(arg(1)),
    repr_(*this, expressions_),
    str_("var" + expressions_[1]->str() + "(" + expressions_[0]->str() + ")")
{
}
//-----------------------------------------------------------------------------
Variable::~Variable()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Variable::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0 = expressions_[0]->operands(name);
  std::vector<Class const*> expr1 = expressions_[1]->operands(name);
  expr0.insert(expr0.end(), expr1.begin(), expr1.end());
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Variable::level_operands(
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
Variable const* Variable::create(Object::repr_t const& repr) const
{
  return new Variable(repr);
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Variable::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
ValueArray const Variable::shape() const
{
  return expressions_[0]->shape();
}

//-----------------------------------------------------------------------------
tuple<Index> const Variable::free_indices() const
{
  return expressions_[0]->free_indices();
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const Variable::index_dimensions() const
{
  return expressions_[0]->index_dimensions();
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const Variable::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Cell const Variable::cell() const
{
  return expressions_[0]->cell();
}

//-----------------------------------------------------------------------------
bool Variable::is_cellwise_constant() const
{
  return expressions_[0]->is_cellwise_constant();
}

//-----------------------------------------------------------------------------
Object::repr_t const& Variable::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Variable::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Variable::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Variable::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  expr.push_back(Expression::create(reprs[0]));
  expr.push_back(new Label(reprs[1]));
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const Variable::fill_expressions(
    Expression const& e, Label const& l)
{
  std::vector<Expression const *> expr;
  expr.push_back(&e);
  expr.push_back(&l);
  return expr;
}
}
