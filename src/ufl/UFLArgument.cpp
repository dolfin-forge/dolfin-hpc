// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLArgument.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Argument::Argument(FiniteElementBase const& fe, dolfin::uint const& c) :
    Expression("Argument"),
    finite_element_(fe),
    count_(c),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "v_" + count_.str() : "v_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Argument::Argument(repr_t const& repr) :
    Expression("Argument", repr),
    finite_element_(*FiniteElementBase::create(arg(0))),
    count_(arg(1)),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "v_" + count_.str() : "v_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Argument::~Argument()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Class const*> const Argument::operands(std::string const& name) const
  {
    std::vector<Class const*> expr0;
    if(name == this->name())
      expr0.push_back(this);
    return expr0;
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<Class const*> > const Argument::level_operands(
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
  Argument const* Argument::create(Object::repr_t const& repr)
  {
    return new Argument(repr);
  }

//-----------------------------------------------------------------------------
  std::vector<Expression const *> const Argument::operands() const
  {
    error("Not yet implemented.");
    std::vector<Expression const*>  expressions_;
    return expressions_;
  }

//-----------------------------------------------------------------------------
  FiniteElementBase const& Argument::element() const
  {
    return finite_element_;
  }

//-----------------------------------------------------------------------------
  type<dolfin::uint> const& Argument::count() const
  {
    return count_;
  }

//-----------------------------------------------------------------------------
  ValueArray const Argument::shape() const
  {
    return finite_element_.value_shape();
  }

//-----------------------------------------------------------------------------
  tuple<Index> const Argument::free_indices() const
  {
    return tuple<Index>();
  }

//-----------------------------------------------------------------------------
  dict<IndexBase, type<dolfin::uint> > const Argument::index_dimensions() const
  {
    return dict<IndexBase, type<dolfin::uint> >();
  }

//-----------------------------------------------------------------------------
  std::vector<std::vector<std::vector<dolfin::real> > > const Argument::evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      ufc::cell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const
  {
    std::cout << "Argument::evaluate " << n << std::endl;
    ufc::finite_element* fe = new dolfin::FiniteElement(element());

    std::cout << "value_rank " << fe->value_rank() << std::endl;
    std::cout << "space_dimension " << fe->space_dimension() << std::endl;
    std::cout << "topological_dimension " << fe->topological_dimension() << std::endl;
    std::cout << "geometric_dimension " << fe->geometric_dimension() << std::endl;

    std::vector<std::vector<std::vector<dolfin::real> > > grads(fe->space_dimension());
    std::vector<std::vector<std::vector<dolfin::real> > > ref_grads(fe->space_dimension());
    for(dolfin::uint s=0; s<fe->space_dimension(); ++s)
    {
      std::cout << "s=" << s << std::endl;
      std::cout << "value_dimension " << fe->value_dimension(0) << std::endl;
      grads[s].resize(q_points.size());
      ref_grads[s].resize(q_points.size());
      for(dolfin::uint qp=0; qp<q_points.size(); ++qp)
      {
        std::cout << "qp=" << qp << ": " << q_points[qp][0] <<  std::endl;
        if(n>0)
        {
          dolfin::uint num_derivatives = fe->value_dimension(0);
          for(dolfin::uint r=0; r<n; ++r)
            num_derivatives *= fe->geometric_dimension();
          grads[s][qp].resize(num_derivatives);
          ref_grads[s][qp].resize(num_derivatives);
          for(dolfin::uint i=0; i<grads[s][qp].size(); ++i)
          {
            grads[s][qp][i] = 0.;
            ref_grads[s][qp][i] = 0.;
          }
          fe->evaluate_basis_derivatives(s, n,
              &grads[s][qp][0], q_points[qp], ref_cell);
          fe->evaluate_reference_basis_derivatives(s, n,
              &ref_grads[s][qp][0], q_points[qp], ref_cell);

          for(dolfin::uint i=0; i<num_derivatives; ++i)
          {
            std::cout << "i= " << i << " " << grads[s][qp][i] << std::endl;
            std::cout << "ref  " << ref_grads[s][qp][i] << std::endl;
          }
        }
        else
        {
          dolfin::uint num_comps = fe->value_dimension(0);
          grads[s][qp].resize(num_comps);
          for(dolfin::uint i=0; i<grads[s][qp].size(); ++i)
          {
            grads[s][qp][i] = 0.;
          }
          fe->evaluate_basis(s, &grads[s][qp][0], q_points[qp], ref_cell);
          for(dolfin::uint i=0; i<num_comps; ++i)
          {
            std::cout << "i= " << i << " " << grads[s][qp][i] << std::endl;
          }
        }
      }
    }
    return grads;
  }

//-----------------------------------------------------------------------------
  Cell const Argument::cell() const
  {
    return finite_element_.cell();
  }

//-----------------------------------------------------------------------------
  bool const Argument::is_cellwise_constant() const
  {
    return false;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& Argument::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Argument::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Argument::display() const
  {
  }

}
