// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/ufl/UFLIndex.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      dolfin::uint const& count) :
    Expression(name),
    count_(count),
    name_(name)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      IndexBase const& index) :
    Expression(name),
    count_(index.count()),
    name_(name)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      repr_t const & repr) :
    Expression(name, repr),
    count_(arg(0)),
    name_(name)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::~IndexBase()
  {
  }

//-----------------------------------------------------------------------------
  type<dolfin::uint> const& IndexBase::count() const
  {
    return count_;
  }

//-----------------------------------------------------------------------------
  void IndexBase::display() const
  {
  }

//-----------------------------------------------------------------------------
//  std::vector<Class const*> const IndexBase::operands(std::string const& name) const
//  {
//    error("Unknown type of ufl::IndexBase: '" + name_ + "'");
//
//    std::vector<Class const*> expr0;
//    return expr0;
//  }

//-----------------------------------------------------------------------------
//  std::vector<std::vector<Class const*> > const IndexBase::level_operands(
//      std::vector<std::vector<Class const*> > const& operands) const
//  {
//    error("Unknown type of ufl::IndexBase: '" + name_ + "'");
//
//    std::vector<std::vector<Class const*> > new_operands0;
//    return new_operands0;
//  }

//-----------------------------------------------------------------------------
  IndexBase const* IndexBase::create(Object::repr_t const& repr)
  {
    std::string name = Class::make_name(repr);
    if (name == "Index")
    {
      return new Index(repr);
    }
    else if (name == "FixedIndex")
    {
      return new FixedIndex(repr);
    }
    else
    {
      error("Unknown type of ufl::IndexBase: '" + name + "'");
    }

    return NULL;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const IndexBase::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const IndexBase::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const IndexBase::free_indices() const
  {
    return tuple<Index>();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const IndexBase::index_dimensions() const
  {
    return dict<IndexBase, type<dolfin::uint> >();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const IndexBase::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Index::Index(type<dolfin::uint> const& count) :
    IndexBase("Index", count),
    repr_(*this, count_),
    str_(count_<10 ? "i_" + count_.str() : "i_{" + count_.str() + "}")
  {
  }

//-----------------------------------------------------------------------------
  Index::Index(repr_t const& repr) :
    IndexBase("Index", repr),
    repr_(*this, count_),
    str_(count_<10 ? "i_" + count_.str() : "i_{" + count_.str() + "}")
  {
  }

//-----------------------------------------------------------------------------
  Index::~Index()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Index::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0;
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Index::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0;
    std::vector<Class const*> obj0;
    obj0.push_back(this);
    new_operands0.push_back(obj0);

    for(dolfin::uint i=0; i<operands.size(); ++i)
      new_operands0.push_back(operands[i]);

    return new_operands0;
  }

//-----------------------------------------------------------------------------
  Index const* Index::create(Object::repr_t const& repr)
  {
    return new Index(repr);
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& Index::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Index::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Index::display() const
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::FixedIndex(dolfin::uint const& count) :
    IndexBase("FixedIndex", count),
    repr_(*this, count_),
    str_(count_.str())
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::FixedIndex(repr_t const& repr) :
    IndexBase("FixedIndex", repr),
    repr_(*this, count_),
    str_(count_.str())
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::~FixedIndex()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const FixedIndex::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0;
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const FixedIndex::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0;
    std::vector<Class const*> obj0;
    obj0.push_back(this);
    new_operands0.push_back(obj0);

    for(dolfin::uint i=0; i<operands.size(); ++i)
      new_operands0.push_back(operands[i]);

    return new_operands0;
  }

//-----------------------------------------------------------------------------
  FixedIndex const* FixedIndex::create(Object::repr_t const& repr)
  {
    return new FixedIndex(repr);
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& FixedIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& FixedIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void FixedIndex::display() const
  {
  }

//-----------------------------------------------------------------------------
//MultiIndex::MultiIndex(tuple<type<dolfin::uint> > const& indices,
//    dict<type<dolfin::uint>, type<dolfin::uint> > const& index_dimensions) :
//    Expression("MultiIndex"),
//    indices_(indices_),
//    index_dimensions_(index_dimensions),
//    repr_(*this, indices_, index_dimensions_),
//    str_("")
//  {
//  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(tuple<IndexBase> const& indices,
      dict<IndexBase, type<dolfin::uint> > const& index_dimensions) :
    Expression("MultiIndex"),
    index_dimensions_(index_dimensions),
    indices_(indices),
//    fixed_indices_(),
    repr_(*this, indices_, index_dimensions_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
//  MultiIndex::MultiIndex(tuple<FixedIndex> const& indices) :
//    Expression("MultiIndex"),
//    index_dimensions_(),
//    indices_(),
//    fixed_indices_(indices),
//    repr_(*this, fixed_indices_, index_dimensions_),
//    str_("")
//  {
//  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(repr_t const& repr) :
    Expression("MultiIndex", repr),
    index_dimensions_(arg(1)),
    indices_(fill_indices(arg(0))),
//    fixed_indices_(fill_fixed_indices(arg(0))),
    repr_(*this, indices_, index_dimensions_),
    str_("")
  {
//    std::stringstream ssrepr;
//    ssrepr << "MultiIndex(";
//    if(indices_.size()>0)
//      ssrepr << indices_.repr() << ", " << index_dimensions_.repr() << ")";
//    else
//      ssrepr << fixed_indices_.repr() << ", " << index_dimensions_.repr() << ")";
//
//    repr_ = ssrepr.str();
  }

//-----------------------------------------------------------------------------
 MultiIndex::~MultiIndex()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const MultiIndex::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = indices_.operands(name);
//    std::vector<Class const*> expr1 = fixed_indices_.operands(name);
    std::vector<Class const*> expr2 = index_dimensions_.operands(name);
//    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    expr0.insert(expr0.end(), expr2.begin(), expr2.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const MultiIndex::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<Class const*> > new_operands0;
      if(indices_.size()>0)
        new_operands0 = indices_.level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands1;
//      if(fixed_indices_.size()>0)
//        new_operands1 = fixed_indices_.level_operands(operands);
    std::vector<std::vector<Class const*> > new_operands2
      = index_dimensions_.level_operands(operands);

    const dolfin::uint size = std::max(new_operands0.size(),
        std::max(new_operands1.size(), new_operands2.size()));
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
          std::cout << new_operands0[i][j]->name() << std::endl;
        }
//      if(i<new_operands1.size())
//        for(dolfin::uint j=0; j<new_operands1[i].size(); ++j)
//        {
//          tmp[i+1].push_back(new_operands1[i][j]);
//          std::cout << new_operands1[i][j]->name() << std::endl;
//        }
      if(i<new_operands2.size())
        for(dolfin::uint j=0; j<new_operands2[i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands2[i][j]);
          std::cout << new_operands2[i][j]->name() << std::endl;
        }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
 MultiIndex const* MultiIndex::create(Object::repr_t const& repr)
{
    return new MultiIndex(repr);
}

//-----------------------------------------------------------------------------
  ValueArray const MultiIndex::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const MultiIndex::free_indices() const
  {
    std::vector<Index const *> indices;
    std::vector<IndexBase const *> base_indices = indices_.operands();
    for(dolfin::uint i=0; i<base_indices.size(); ++i)
    {
      if(base_indices[i]->name() == "Index")
        indices.push_back(dynamic_cast<Index const *>(base_indices[i]));
    }

    return indices;
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const MultiIndex::index_dimensions() const
  {
    return index_dimensions_;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const MultiIndex::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
 Object::repr_t const& MultiIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
 std::string const& MultiIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
 void MultiIndex::display() const
  {
  }

//-----------------------------------------------------------------------------
  tuple<IndexBase> const MultiIndex::fill_indices(repr_t const& repr)
  {
    if(index_dimensions_.size() > 0)
      return tuple<IndexBase>(repr);

    return tuple<IndexBase>();
  }

//-----------------------------------------------------------------------------
//  tuple<FixedIndex> const MultiIndex::fill_fixed_indices(repr_t const& repr)
//  {
//    if(index_dimensions_.size() == 0)
//      return tuple<FixedIndex>(repr);
//
//    return tuple<FixedIndex>();
//  }

//-----------------------------------------------------------------------------
 std::vector<Expression const *> const MultiIndex::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

}
