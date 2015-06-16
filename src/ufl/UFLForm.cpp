// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/ufl/UFLForm.h>

#include <dolfin/common/types.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Form::Form(List const& list) :
    Class("Form"),
    list_(list),
    form_data_(new FormData(*this)),
    repr_(*this, list_),
    str_("\n  +  " + list_.str())
//    is_preprocessed_(false)
  {
    std::cout << "in Form" << std::endl;
    std::cout << repr_ << std::endl;
    std::cout << form_data_->n_cell_elements() << " "
      <<  form_data_->n_extf_elements() << " "
      <<  form_data_->n_intf_elements() << std::endl;
  }

//-----------------------------------------------------------------------------
  Form::Form(repr_t const& repr) :
    Class("Form", repr),
    list_(arg(0)),
    form_data_(new FormData(*this)),
    repr_(*this, list_),
    str_("\n  +  " + list_.str())
//    is_preprocessed_(false)
  {
    std::cout << "in Form" << std::endl;
    std::cout << repr_ << std::endl;
    std::cout << form_data_->n_cell_elements() << " "
      <<  form_data_->n_extf_elements() << " "
      <<  form_data_->n_intf_elements() << std::endl;
  }

//-----------------------------------------------------------------------------
  Form::Form(ufc::form const& ufc_form) :
    Class("Form", ufc_form.signature()),
    list_(arg(0)),
    form_data_(new FormData(*this)),
    repr_(*this, list_),
    str_("\n  +  " + list_.str())
//    is_preprocessed_(false)
  {
    std::cout << "in Form" << std::endl;
    std::cout << repr_ << std::endl;
    std::cout << form_data_->n_cell_elements() << " "
      <<  form_data_->n_extf_elements() << " "
      <<  form_data_->n_intf_elements() << std::endl;
  }

//-----------------------------------------------------------------------------
  Form::~Form()
  {
  }

//-----------------------------------------------------------------------------
  Cell const& Form::cell() const //Does this function make sense?
  {
    std::vector<Integral const *> ints = cell_integrals().operands();
    for(uint i = 0; i < ints.size(); ++i)
    {
      std::vector<Expression const *> integrands = ints[i]->integrand();
      for(uint j = 0; j < integrands.size(); ++j)
        return integrands[j]->cell();
    }

    error("No cell found in an integrand.");
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::integrals(
      MeasureDomain::Type const& measure_type) const
  {
    //get vector of integrals
    std::vector<Integral const *> const& integrals = list_.get_integrals();
    std::vector<Integral const *> integrals_to_return;

    for(dolfin::uint i = 0; i<integrals.size(); ++i)
    {
      if(integrals[i]->measure().measure_type() == measure_type)
        integrals_to_return.push_back(integrals[i]);
    }

    return tuple<Integral> (integrals_to_return);
  }

//-----------------------------------------------------------------------------
  tuple<Measure> const Form::measures(
      MeasureDomain::Type const& measure_type) const
  {
    //get vector of integrals
    std::vector<Integral const *> const& integrals = list_.get_integrals();
    std::vector<Measure const *> measures_to_return;

    for(dolfin::uint i = 0; i<integrals.size(); ++i)
      if(integrals[i]->measure().measure_type() == measure_type)
        measures_to_return.push_back(&(integrals[i]->measure()));

    return tuple<Measure> (measures_to_return);
  }

//-----------------------------------------------------------------------------
  std::vector<std::pair<MeasureDomain::Type, dolfin::uint> > const Form::domains(
      MeasureDomain::Type const& measure_type) const
  {
    //get vector of integrals
    std::vector<Integral const *> const& integrals = list_.get_integrals();
    std::vector<std::pair<MeasureDomain::Type, dolfin::uint> > domains;

    for(dolfin::uint i = 0; i<integrals.size(); ++i)
      if(integrals[i]->measure().measure_type() == measure_type)
        domains.push_back(std::make_pair(measure_type, integrals[i]->measure().measure_domain_id()));

    return domains;
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::cell_integrals() const
  {
    return integrals(MeasureDomain::cell);
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::exterior_facet_integrals() const
  {
    return integrals(MeasureDomain::exterior_facet);
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::interior_facet_integrals() const
  {
    return integrals(MeasureDomain::interior_facet);
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::macro_cell_integrals() const
  {
    return integrals(MeasureDomain::macro_cell);
  }

//-----------------------------------------------------------------------------
  tuple<Integral> const Form::surface_integrals() const
  {
    return integrals(MeasureDomain::surface);
  }

//-----------------------------------------------------------------------------
  FormData const& Form::form_data() const
  {
    return *form_data_;
  }

//-----------------------------------------------------------------------------
//  bool const Form::is_preprocessed() const
//  {
//    return is_preprocessed_;
//  }

//-----------------------------------------------------------------------------
  Object::repr_t const& Form::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const& Form::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Form::display() const
  {
  }
}
