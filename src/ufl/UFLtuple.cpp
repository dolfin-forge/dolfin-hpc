// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLTuple.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Tuple::Tuple(tuple<Expression> const& t) :
    Expression("Tuple"),
    t_(t)
  {
    // Create string representation
    std::stringstream ssrepr;
    std::stringstream ssstr;
    ssrepr << "Tuple(*" << t_.repr() << ")";
    ssstr << "Tuple(*" << t_.str() << ")";

    repr_ = ssrepr.str();
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  Tuple::Tuple(repr_t const & repr) :
    Expression("Tuple", repr),
    t_(arg(0))
  {
    // Create string representation
    std::stringstream ssrepr;
    std::stringstream ssstr;
    ssrepr << "Tuple(*" << t_.repr() << ")";
    ssstr << "Tuple(*" << t_.str() << ")";

    repr_ = ssrepr.str();
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  Tuple::~Tuple()
  {
  }
  
//-----------------------------------------------------------------------------
  std::vector<Class const*> const Tuple::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0;
    for(dolfin::uint i=0; i<operands().size(); ++i)
    {
     std::vector<Class const*> expr1 = operands()[i]->operands(name);
     expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    }
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Tuple::level_operands(
      std::vector<std::vector<Class const*> > const& ops) const
  {
    std::vector<std::vector<std::vector<Class const*> > > new_operands(operands().size());

    dolfin::uint size = 0;
    for(dolfin::uint i=0; i<operands().size(); ++i)
    {
      new_operands[i] = operands()[i]->level_operands(ops);
      size = std::max(size, (dolfin::uint)new_operands[i].size());
    }

    std::vector<std::vector<Class const*> > tmp(size+1);
    std::vector<Class const*> obj0;
    obj0.push_back(this);

    tmp[0] = obj0;
    for(dolfin::uint i=0; i<tmp.size()-1; ++i)
    {
      for(dolfin::uint k=0; k<new_operands.size(); ++k)
      {
        if(i<new_operands[k].size())
          for(dolfin::uint j=0; j<new_operands[k][i].size(); ++j)
          {
            tmp[i+1].push_back(new_operands[k][i][j]);
//            std::cout << new_operands[k][i][j]->name() << std::endl;
          }
      }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Tuple::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const Tuple::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Tuple::display() const
  {
  }

//-----------------------------------------------------------------------------
  Tuple const* Tuple::create(Object::repr_t const& repr) const
  {
      return new Tuple(repr);
  }
  
//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Tuple::operands() const
  {
    return t_.operands();  
  }

//-----------------------------------------------------------------------------
  ValueArray const Tuple::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Tuple::free_indices() const
  {
    return tuple<Index>();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Tuple::index_dimensions() const
  {
    return dict<IndexBase, type<dolfin::uint> >();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Tuple::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell, 
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    for(dolfin::uint i=0; i<operands().size(); ++i)
      operands()[i]->evaluate(n, tensor, ref_cell, q_points, coordinates);
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

}
