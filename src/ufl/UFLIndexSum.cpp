// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLIndexSum.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
IndexSum::IndexSum(Expression const& summand, MultiIndex const& index) :
    Expression("IndexSum"),
    expressions_(fill_expressions(summand, index)),
    repr_(*this, expressions_),
    str_("sum_{" + index.str() + "} " + summand.str() + " ")
{
}

//-----------------------------------------------------------------------------
IndexSum::IndexSum(Expression const& expression,
                   tuple<IndexBase> const& indices,
                   dict<IndexBase, type<dolfin::uint> > const& index_dimensions) :
    Expression("IndexSum"),
    expressions_(fill_expressions(expression, indices, index_dimensions)),
    repr_(*this, expressions_),
    str_("sum_{" + expressions_[1]->str() + "} " + expressions_[0]->str() + " ")
{
}

//-----------------------------------------------------------------------------
IndexSum::IndexSum(repr_t const & repr) :
    Expression("IndexSum", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("sum_{" + expressions_[1]->str() + "} " + expressions_[0]->str() + " ")
{
}

//-----------------------------------------------------------------------------
IndexSum::~IndexSum()
= default;

//-----------------------------------------------------------------------------
std::vector<Class const*> const IndexSum::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0 = expressions_[0]->operands(name);
  std::vector<Class const*> expr1 = expressions_[1]->operands(name);
  expr0.insert(expr0.end(), expr1.begin(), expr1.end());
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const IndexSum::level_operands(
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
IndexSum const* IndexSum::create(Object::repr_t const& repr)
{
  return new IndexSum(repr);
}

//-----------------------------------------------------------------------------
MultiIndex const * IndexSum::index() const
{
  return new MultiIndex(expressions_[1]->repr());
}

//-----------------------------------------------------------------------------
dolfin::uint IndexSum::dimension() const
{
  std::vector<Index const *> const& f_indices =
      expressions_[1]->free_indices().operands();

  if (f_indices.size() != 1) error("Expecting a single Index only.");

  Index const& ind = *f_indices[0];
  std::vector<std::pair<IndexBase const *, type<dolfin::uint> const *> > const& ind_dim =
      index_dimensions().map();

  for (dolfin::uint i = 0; i < ind_dim.size(); ++i)
    if (*ind_dim[i].first == ind) return *ind_dim[i].second;

  return -1;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const IndexSum::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
ValueArray const IndexSum::shape() const
{
  return expressions_[0]->shape();
}

//-----------------------------------------------------------------------------
tuple<Index> const IndexSum::free_indices() const
{
  std::vector<Index const *> const& sum_indices =
      expressions_[0]->free_indices().operands();

  std::vector<Index const *> const& f_indices =
      expressions_[1]->free_indices().operands();

  if (f_indices.size() != 1) error("Expecting a single Index only.");

  std::vector<Index const *> new_indices;

  for (dolfin::uint i = 0; i < sum_indices.size(); ++i)
    if (sum_indices[i] != f_indices[0]) new_indices.push_back(sum_indices[i]);

  return tuple<Index>(new_indices);
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const IndexSum::index_dimensions() const
{
  return expressions_[0]->index_dimensions();
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const IndexSum::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    ufc::cell const& ref_cell, std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  std::cout << "IndexSum::evaluate" << n << std::endl;
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0 =
      expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
  dolfin::uint tensor_rank = expressions_[1]->index_dimensions().size();

  std::vector<std::vector<std::vector<dolfin::real> > > result(
      new_vals0.size());
  for (dolfin::uint i = 0; i < new_vals0.size(); ++i)
  {
    result[i].resize(new_vals0[i].size());
    for (dolfin::uint j = 0; j < new_vals0[i].size(); ++j)
    {
      result[i][j].resize(tensor_rank);
      for (dolfin::uint k = 0; k < new_vals0[i][j].size(); ++k)
        result[i][j][k] = new_vals0[i][j][k];  //TODO
    }
  }

  for (dolfin::uint i = 0; i < new_vals0.size(); ++i)
    for (dolfin::uint j = 0; j < new_vals0[i].size(); ++j)
      for (dolfin::uint k = 0; k < new_vals0[i][j].size(); ++k)
        std::cout << "(" << i << "," << j << "," << k << ") = "
                  << new_vals0[i][j][k] << std::endl;
  std::cout << std::endl;
//    for(dolfin::uint i=0; i<new_vals1.size(); ++i)
//      for(dolfin::uint j=0; j<new_vals1[i].size(); ++j)
//        for(dolfin::uint k=0; k<new_vals1[i][j].size(); ++k)
//          std::cout << "(" << i << "," << j << "," << k << ") = " << new_vals1[i][j][k] << std::endl;
  std::cout << "IndexSum::evaluate END" << n << std::endl;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Object::repr_t const& IndexSum::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& IndexSum::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void IndexSum::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const IndexSum::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  expr.push_back(Expression::create(reprs[0]));
  expr.push_back(new MultiIndex(reprs[1]));
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const IndexSum::fill_expressions(
    Expression const& expression, MultiIndex const& multi_index)
{
  std::vector<Expression const *> expr;
  expr.push_back(&expression);
  expr.push_back(&multi_index);
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const IndexSum::fill_expressions(
    Expression const& expression, tuple<IndexBase> const& indices,
    dict<IndexBase, type<dolfin::uint> > const& index_dimensions)
{
  std::vector<Expression const *> expr;
  expr.push_back(&expression);
  expr.push_back(new MultiIndex(indices, index_dimensions));
  return expr;
}

}
