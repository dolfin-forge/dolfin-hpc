// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIntegral.h>
#include <dolfin/ufl/UFLrepr.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
MeasureDomain::MappingReprToType const MeasureDomain::__init_mapping_repr_to_type()
{
  MappingReprToType m;
  m.insert(MappingReprToTypeItem(Object::repr_t("None"), MeasureDomain::None));
  m.insert(MappingReprToTypeItem(Object::repr_t("'cell'"), MeasureDomain::cell));
  m.insert(MappingReprToTypeItem(Object::repr_t("'exterior_facet'"), MeasureDomain::exterior_facet));
  m.insert(MappingReprToTypeItem(Object::repr_t("'interior_facet'"), MeasureDomain::interior_facet));
  m.insert(MappingReprToTypeItem(Object::repr_t("'macro_cell'"), MeasureDomain::macro_cell));
  m.insert(MappingReprToTypeItem(Object::repr_t("'surface'"), MeasureDomain::surface));
  return m;
}

//-----------------------------------------------------------------------------
MeasureDomain::MappingTypeToRepr const MeasureDomain::__init_mapping_type_to_repr()
{
  MappingTypeToRepr m;
  m.insert(MappingTypeToReprItem(MeasureDomain::None, Object::repr_t("None")));
  m.insert(MappingTypeToReprItem(MeasureDomain::cell, Object::repr_t("'cell'")));
  m.insert(MappingTypeToReprItem(MeasureDomain::exterior_facet, Object::repr_t("'exterior_facet'")));
  m.insert(MappingTypeToReprItem(MeasureDomain::interior_facet, Object::repr_t("'interior_facet'")));
  m.insert(MappingTypeToReprItem(MeasureDomain::macro_cell, Object::repr_t("'macro_cell'")));
  m.insert(MappingTypeToReprItem(MeasureDomain::surface, Object::repr_t("'surface'")));
  return m;
}

//-----------------------------------------------------------------------------
MeasureDomain::MeasureDomain(MeasureDomain::Type const& t) :
    ufl::type<std::string>(MeasureDomain::type_repr(t)),
    type_(t)
{
}

//-----------------------------------------------------------------------------
MeasureDomain::MeasureDomain(repr_t const& repr) :
    ufl::type<std::string>(repr),
    type_(MeasureDomain::repr_type(repr))
{
}

//-----------------------------------------------------------------------------
MeasureDomain::~MeasureDomain()
{
}

//-----------------------------------------------------------------------------
std::string const MeasureDomain::type_repr(MeasureDomain::Type const& t)
{
  return MappingTypeToRepr().find(t)->second;
}

//-----------------------------------------------------------------------------
MeasureDomain::Type const MeasureDomain::repr_type(repr_t const& repr)
{
  return MappingReprToType().find(repr)->second;
}

//-----------------------------------------------------------------------------
MeasureDomain::Type const MeasureDomain::type() const
{
  return type_;
}

//-----------------------------------------------------------------------------
void MeasureDomain::display() const
{
//  ufl::type<std::string>::display();
//  std::cout << std::setw(24) << "dimension" << " = " << this->dim()
//      << std::endl;
//  std::cout << std::setw(24) << "facet" << " = "
//      << Domain(this->facet()).str() << std::endl;
//  std::cout << std::setw(24) << "num_facets" << " = " << this->num_facets()
//      << std::endl;
//  std::cout << std::endl;
}

//-----------------------------------------------------------------------------
  Measure::Measure(MeasureDomain::Type const& measure_type, 
      dolfin::uint measure_id) :
    Class("Measure"),
    measure_(measure_type),
    measure_id_(measure_id)
  {
  }

//-----------------------------------------------------------------------------
  Measure::Measure(repr_t const & repr) :
    Class("Measure", repr),
    measure_(),
    measure_id_(0)
  {
  }

//-----------------------------------------------------------------------------
  Measure::~Measure()
  {
  }

//-----------------------------------------------------------------------------
  MeasureDomain::Type const& Measure::measure_type() const
  {
    return measure_;
  }

//-----------------------------------------------------------------------------
  dolfin::uint Measure::measure_id() const
  {
    return measure_id_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Measure::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Measure::str() const
  {
    return str_;
  }
  
//-----------------------------------------------------------------------------
  void Measure::display() const
  {
  }


//-----------------------------------------------------------------------------
  Integral::Integral(Expression const & integrand, Measure const & measure) :
    Class("Integral"),
    integrand_(integrand),
    measure_(measure) 
  {
  }

//-----------------------------------------------------------------------------
  Integral::Integral(repr_t const & repr) :
    Class("Integral", repr),
    integrand_(arg(0)),
    measure_(arg(1))
  {
  }

//-----------------------------------------------------------------------------
  Integral::~Integral()
  {
  }

//-----------------------------------------------------------------------------
  Expression const& Integral::integrand() const
  {
    return integrand_;
  }

//-----------------------------------------------------------------------------
  Measure const& Integral::measure() const
  {
    return measure_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Integral::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Integral::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Integral::display() const
  {
  }


}
