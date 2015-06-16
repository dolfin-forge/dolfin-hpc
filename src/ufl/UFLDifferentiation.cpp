// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/ufl/UFLDifferentiation.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  CoefficientDerivative::CoefficientDerivative(type<dolfin::real> const& integrand,
      Tuple const& coefficients, Tuple const& arguments,
      dict<Data, type<dolfin::uint> > const& coeff_derivatives) :
    Expression("CoefficientDerivative"),
    integrand_(integrand),
    expressions_(fill_expressions(coefficients, arguments)),
    coeff_derivatives_(coeff_derivatives),
    repr_(*this, integrand_, coefficients, arguments, coeff_derivatives_),
    str_("d/dfj { " + integrand_.str() + " } with fh=" + coefficients.str() + ", dfh/dfj = " + arguments.str() + ". and coefficient derivatives " + coeff_derivatives_.str())
  {
  }

//-----------------------------------------------------------------------------
  CoefficientDerivative::CoefficientDerivative(type<dolfin::real> const& integrand,
      Tuple const& coefficients, Tuple const& arguments,
      Data const& coeff_derivatives) :
    Expression("CoefficientDerivative"),
    integrand_(integrand),
    expressions_(fill_expressions(coefficients, arguments)),
    coeff_derivatives_(dict<Data, type<dolfin::uint> > (std::pair<Data const*, type<dolfin::uint> const *>(&coeff_derivatives, new type<dolfin::uint>(1)))),
    repr_(*this, integrand_, coefficients, arguments, coeff_derivatives_),
    str_("d/dfj { " + integrand_.str() + " } with fh=" + coefficients.str() + ", dfh/dfj = " + arguments.str() + ". and coefficient derivatives " + coeff_derivatives_.str())
  {
  }

//-----------------------------------------------------------------------------
  CoefficientDerivative::CoefficientDerivative(repr_t const & repr) :
    Expression("CoefficientDerivative", repr),
    integrand_(arg(0)),
    expressions_(fill_expressions(args())),
    coeff_derivatives_(arg(3)),
    repr_(*this, integrand_, *expressions_[0], *expressions_[1], coeff_derivatives_),
    str_("d/dfj { " + integrand_.str() + " } with fh=" + expressions_[0]->str() + ", dfh/dfj = "
        + expressions_[1]->str() + ". and coefficient derivatives " + coeff_derivatives_.str())
  {
  }

//-----------------------------------------------------------------------------
  CoefficientDerivative::~CoefficientDerivative()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const CoefficientDerivative::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const CoefficientDerivative::level_operands(
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
  CoefficientDerivative const* CoefficientDerivative::create(Object::repr_t const& repr)
  {
    return new CoefficientDerivative(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const CoefficientDerivative::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const CoefficientDerivative::shape() const
  {
    return ValueArray();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const CoefficientDerivative::free_indices() const
  {
    return tuple<Index>();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const CoefficientDerivative::index_dimensions() const
  {
    return dict<IndexBase, type<dolfin::uint> >();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const CoefficientDerivative::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    //FIXME implement this
    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& CoefficientDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& CoefficientDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void CoefficientDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const CoefficientDerivative::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(new Tuple(reprs[1]));
    expr.push_back(new Tuple(reprs[2]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const CoefficientDerivative::fill_expressions(Tuple const& t1, Tuple const& t2)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&t1);
    expr.push_back(&t2);
    return expr;
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::SpatialDerivative(Expression const& expression, IndexBase const& index,
      type<dolfin::uint> const& index_dimension) :
    Expression("SpatialDerivative"),
    expressions_(fill_expressions(expression,
          MultiIndex(tuple<IndexBase>(index), dict<IndexBase, type<dolfin::uint> > (index, index_dimension)))),
    repr_(*this, expressions_),
    str_("d/dx_" + index.str() + " " + expression.str())
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::SpatialDerivative(Expression const& expression, MultiIndex const& index) :
    Expression("SpatialDerivative"),
    expressions_(fill_expressions(expression, index)),
    repr_(*this, expressions_),
    str_("d/dx_" + index.str() + " " + expression.str())
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::SpatialDerivative(repr_t const& repr) :
    Expression("SpatialDerivative", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("d/dx_" + expressions_[1]->str() + " " + expressions_[0]->str())
  {
  }

//-----------------------------------------------------------------------------
  SpatialDerivative::~SpatialDerivative()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const SpatialDerivative::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const SpatialDerivative::level_operands(
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
  SpatialDerivative const* SpatialDerivative::create(Object::repr_t const& repr)
  {
    return new SpatialDerivative(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const SpatialDerivative::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const SpatialDerivative::shape() const
  {
    return expressions_[0]->shape();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const SpatialDerivative::free_indices() const
  {
    return expressions_[1]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const SpatialDerivative::index_dimensions() const
  {
    return expressions_[1]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const SpatialDerivative::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::cout << "SpatialDerivative::evaluate " << n << std::endl;
    n++;
    std::vector<Index const*> const f_indices =
      expressions_[1]->free_indices().operands();
    std::vector<std::pair<IndexBase const *, type<dolfin::uint> const*> > const index_dim =
      expressions_[1]->index_dimensions().map();

    for(dolfin::uint i=0; i<f_indices.size(); ++i)
    {
      std::cout << "i=" << i << ": " << f_indices[i]->count();
      std::cout << "  " << *index_dim[i].second << std::endl;
    }

    std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0
      = expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);

    for(dolfin::uint i=0; i<new_vals0.size(); ++i)
      for(dolfin::uint j=0; j<new_vals0[i].size(); ++j)
        for(dolfin::uint k=0; k<new_vals0[i][j].size(); ++k)
          std::cout << "(" << i << "," << j << "," << k << ") = " << new_vals0[i][j][k] << std::endl;

    std::cout << "SpatialDerivative::evaluate END" << n << std::endl;
    return new_vals0;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& SpatialDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& SpatialDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void SpatialDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const SpatialDerivative::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(new MultiIndex(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const SpatialDerivative::fill_expressions(Expression const& e, MultiIndex const& i)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    expr.push_back(&i);
    return expr;
  }

//-----------------------------------------------------------------------------
  VariableDerivative::VariableDerivative(Expression const& expression, Variable const& variable) :
    Expression("VariableDerivative"),
    expressions_(fill_expressions(expression, variable)),
    repr_(*this, expressions_),
    str_("d/d[" + variable.str() + "] " + expression.str())
  {
  }

//-----------------------------------------------------------------------------
  VariableDerivative::VariableDerivative(repr_t const & repr) :
    Expression("VariableDerivative", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("d/d[" + expressions_[1]->str() + "] " + expressions_[0]->str())
  {
  }

//-----------------------------------------------------------------------------
  VariableDerivative::~VariableDerivative()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const VariableDerivative::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    std::vector<Class const*> expr1 = expressions_[1]->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const VariableDerivative::level_operands(
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
  VariableDerivative const* VariableDerivative::create(Object::repr_t const& repr)
  {
    return new VariableDerivative(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const VariableDerivative::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const VariableDerivative::shape() const
  {
    ValueArray const& shape_0 = expressions_[0]->shape();
    ValueArray const& shape_1 = expressions_[1]->shape();

    ValueArray return_shape = shape_0;
    for(dolfin::uint i=0; i<shape_1.size(); ++i)
      return_shape.push_back(shape_1[i]);
    return return_shape;
  }

//-----------------------------------------------------------------------------
  tuple<Index> const VariableDerivative::free_indices() const
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
  dict<IndexBase, type<dolfin::uint> > const VariableDerivative::index_dimensions() const
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
  std::vector<std::vector<std::vector<dolfin::real> > > const VariableDerivative::evaluate(
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
  Object::repr_t const& VariableDerivative::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& VariableDerivative::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void VariableDerivative::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const VariableDerivative::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    expr.push_back(new Variable(reprs[1]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const VariableDerivative::fill_expressions(Expression const& e, Variable const& v)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    expr.push_back(&v);
    return expr;
  }

//-----------------------------------------------------------------------------
  Grad::Grad(Expression const& expression) :
    Expression("Grad"),
    expressions_(fill_expressions(expression)),
    repr_(*this, expressions_),
    str_("grad(" + expression.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Grad::Grad(repr_t const & repr) :
    Expression("Grad", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("grad(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Grad::~Grad()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Grad::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Grad::level_operands(
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
    return new_operands0;
  }

//-----------------------------------------------------------------------------
  Grad const* Grad::create(Object::repr_t const& repr)
  {
    return new Grad(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Grad::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const Grad::shape() const
  {
    ValueArray const& shape_array = expressions_[0]->shape();
    ValueArray return_array = shape_array;
    return_array.push_back(expressions_[0]->geometric_dimension());

    return return_array;
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Grad::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Grad::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Grad::evaluate(
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
  Object::repr_t const& Grad::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Grad::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Grad::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Grad::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Grad::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

//-----------------------------------------------------------------------------
  Div::Div(Expression const& expression) :
    Expression("Div"),
    expressions_(fill_expressions(expression)),
    repr_(*this, expressions_),
    str_("div(" + expression.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Div::Div(repr_t const & repr) :
    Expression("Div", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("div(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Div::~Div()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Div::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Div::level_operands(
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
  Div const* Div::create(Object::repr_t const& repr)
  {
    return new Div(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Div::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const Div::shape() const
  {
    ValueArray const& shape_array = expressions_[0]->shape();
    ValueArray return_array = shape_array;
    return_array.pop_back();

    return return_array;
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Div::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Div::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Div::evaluate(
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
  Object::repr_t const& Div::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Div::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Div::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Div::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Div::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

//-----------------------------------------------------------------------------
  NablaGrad::NablaGrad(Expression const& expression) :
    Expression("NablaGrad"),
    expressions_(fill_expressions(expression)),
    repr_(*this, expressions_),
    str_("nabla_grad(" + expression.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaGrad::NablaGrad(repr_t const & repr) :
    Expression("NablaGrad", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("nabla_grad(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaGrad::~NablaGrad()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const NablaGrad::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const NablaGrad::level_operands(
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
  NablaGrad const* NablaGrad::create(Object::repr_t const& repr)
  {
    return new NablaGrad(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaGrad::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const NablaGrad::shape() const
  {
    ValueArray const& shape_array = expressions_[0]->shape();
    ValueArray return_array(expressions_[0]->geometric_dimension());
    for(dolfin::uint i=0; i<shape_array.size(); ++i)
      return_array.push_back(shape_array[i]);

    return return_array;
  }

//-----------------------------------------------------------------------------
  tuple<Index> const NablaGrad::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const NablaGrad::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const NablaGrad::evaluate(
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
  Object::repr_t const& NablaGrad::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& NablaGrad::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void NablaGrad::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaGrad::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaGrad::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

//-----------------------------------------------------------------------------
  NablaDiv::NablaDiv(Expression const& expression) :
    Expression("NablaDiv"),
    expressions_(fill_expressions(expression)),
    repr_(*this, expressions_),
    str_("nabla_div(" + expression.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaDiv::NablaDiv(repr_t const & repr) :
    Expression("NablaDiv", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("nabla_div(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  NablaDiv::~NablaDiv()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const NablaDiv::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const NablaDiv::level_operands(
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
  NablaDiv const* NablaDiv::create(Object::repr_t const& repr)
  {
    return new NablaDiv(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaDiv::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const NablaDiv::shape() const
  {
    ValueArray const& shape_array = expressions_[0]->shape();
    ValueArray return_array;

    for(dolfin::uint i=1; i<shape_array.size(); ++i)
      return_array.push_back(shape_array[i]);

    return return_array;
  }

//-----------------------------------------------------------------------------
  tuple<Index> const NablaDiv::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const NablaDiv::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const NablaDiv::evaluate(
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
  Object::repr_t const& NablaDiv::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& NablaDiv::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void NablaDiv::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaDiv::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const NablaDiv::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

//-----------------------------------------------------------------------------
  Curl::Curl(Expression const& expression) :
    Expression("Curl"),
    expressions_(fill_expressions(expression)),
    repr_(*this, expressions_),
    str_("curl(" + expression.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Curl::Curl(repr_t const & repr) :
    Expression("Curl", repr),
    expressions_(fill_expressions(args())),
    repr_(*this, expressions_),
    str_("curl(" + expressions_[0]->str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Curl::~Curl()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Curl::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0 = expressions_[0]->operands(name);
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Curl::level_operands(
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
  Curl const* Curl::create(Object::repr_t const& repr)
  {
    return new Curl(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Curl::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  ValueArray const Curl::shape() const
  {
    return expressions_[0]->shape();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Curl::free_indices() const
  {
    return expressions_[0]->free_indices();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Curl::index_dimensions() const
  {
    return expressions_[0]->index_dimensions();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Curl::evaluate(
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
  Object::repr_t const& Curl::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Curl::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Curl::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Curl::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Expression const *> expr;
    expr.push_back(Expression::create(reprs[0]));
    return expr;
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Curl::fill_expressions(Expression const& e)
  {
    std::vector<Expression const *> expr;
    expr.push_back(&e);
    return expr;
  }

}
