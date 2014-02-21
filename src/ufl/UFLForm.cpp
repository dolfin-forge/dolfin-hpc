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
    repr_(*this, list_),
    str_("\n  +  " + list_.str()),
    is_preprocessed_(false)
  {
  }

//-----------------------------------------------------------------------------
  Form::Form(repr_t const & repr) :
    Class("Form", repr),
    list_(arg(0)),
    repr_(*this, list_),
    str_("\n  +  " + list_.str()),
    is_preprocessed_(false)
  {
    std::cout << "C Form " << repr << std::endl;
  }

//-----------------------------------------------------------------------------
  Form::~Form()
  {
  }
  
//-----------------------------------------------------------------------------
  Cell const& Form::cell() const
  {
//    return integral_.integrand().cell();  
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::integrals(
      MeasureDomain::Type const& measure_type) const
  {
//    if(measure_type == MeasureDomain::None)
//      return integrals_;  

    //tuple(itg for itg in self._integrals if itg.measure().domain_type() == domain_type)
//     return integrals_;
  }

//-----------------------------------------------------------------------------
  std::vector<Measure> const Form::measures(
      MeasureDomain::Type const& measure_type) const
  {
    //tuple(itg.measure() for itg in self.integrals(domain_type))
    std::vector<Measure> measures;
    return measures;
  }

//-----------------------------------------------------------------------------
  std::vector<MeasureDomain::Type> const Form::domains(
      MeasureDomain::Type const& measure_type) const
  {
    //tuple((m.domain_type(), m.domain_id()) for m in self.measures(domain_type))
    std::vector<MeasureDomain::Type> domains;
    return domains;
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::cell_integrals() const
  {
    return integrals(MeasureDomain::cell);
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::exterior_facet_integrals() const
  {
    return integrals(MeasureDomain::exterior_facet);
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::interior_facet_integrals() const
  {
    return integrals(MeasureDomain::interior_facet);
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::macro_cell_integrals() const
  {
    return integrals(MeasureDomain::macro_cell);
  }

//-----------------------------------------------------------------------------
  std::vector<Integral> const& Form::surface_integrals() const
  {
    return integrals(MeasureDomain::surface);
  }

//-----------------------------------------------------------------------------
  bool const Form::is_preprocessed() const
  {
    return is_preprocessed_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Form::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Form::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Form::display() const
  {
  }
}
