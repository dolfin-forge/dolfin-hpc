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
MeasureData::MappingReprToType MeasureData::__init_mapping_repr_to_type()
{
  MappingReprToType m;
  m.insert(MappingReprToTypeItem(Object::repr_t("None"), MeasureData::None));
  return m;
}

//-----------------------------------------------------------------------------
MeasureData::MappingTypeToRepr MeasureData::__init_mapping_type_to_repr()
{
  MappingTypeToRepr m;
  m.insert(MappingTypeToReprItem(MeasureData::None, Object::repr_t("None")));
  return m;
}

//-----------------------------------------------------------------------------
MeasureData::MeasureData(MeasureData::Type const& t) :
    ufl::type<std::string>(MeasureData::type_repr(t)),
    type_(t)
{
}

//-----------------------------------------------------------------------------
MeasureData::MeasureData(repr_t const& repr) :
    ufl::type<std::string>(repr),
    type_(MeasureData::repr_type(repr))
{
}

//-----------------------------------------------------------------------------
MeasureData::~MeasureData()
{
}

//-----------------------------------------------------------------------------
std::string MeasureData::type_repr(MeasureData::Type const& t)
{
  return MappingToRepr().find(t)->second;
}

//-----------------------------------------------------------------------------
MeasureData::Type MeasureData::repr_type(repr_t const& repr)
{
  return MappingToType().find(repr)->second;
}

//-----------------------------------------------------------------------------
MeasureData::Type MeasureData::type() const
{
  return type_;
}

//-----------------------------------------------------------------------------
void MeasureData::display() const
{
}

//-----------------------------------------------------------------------------
MeasureDomain::MappingReprToType MeasureDomain::__init_mapping_repr_to_type()
{
  MappingReprToType m;
  m.insert(MappingReprToTypeItem(Object::repr_t("None"), MeasureDomain::None));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'cell'"), MeasureDomain::cell));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'exterior_facet'"),
                            MeasureDomain::exterior_facet));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'interior_facet'"),
                            MeasureDomain::interior_facet));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'macro_cell'"),
                            MeasureDomain::macro_cell));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'surface'"),
                            MeasureDomain::surface));
  m.insert(MappingReprToTypeItem(Object::repr_t("'dx'"), MeasureDomain::cell));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'ds'"),
                            MeasureDomain::exterior_facet));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'dS'"),
                            MeasureDomain::interior_facet));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'dE'"), MeasureDomain::macro_cell));
  m.insert(
      MappingReprToTypeItem(Object::repr_t("'dc'"), MeasureDomain::surface));
  return m;
}

//-----------------------------------------------------------------------------
MeasureDomain::MappingTypeToRepr MeasureDomain::__init_mapping_type_to_repr()
{
  MappingTypeToRepr m;
  m.insert(MappingTypeToReprItem(MeasureDomain::None, Object::repr_t("None")));
  m.insert(
      MappingTypeToReprItem(MeasureDomain::cell, Object::repr_t("'cell'")));
  m.insert(
      MappingTypeToReprItem(MeasureDomain::exterior_facet,
                            Object::repr_t("'exterior_facet'")));
  m.insert(
      MappingTypeToReprItem(MeasureDomain::interior_facet,
                            Object::repr_t("'interior_facet'")));
  m.insert(
      MappingTypeToReprItem(MeasureDomain::macro_cell,
                            Object::repr_t("'macro_cell'")));
  m.insert(
      MappingTypeToReprItem(MeasureDomain::surface,
                            Object::repr_t("'surface'")));
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
std::string MeasureDomain::type_repr(MeasureDomain::Type const& t)
{
  return MappingToRepr().find(t)->second;
}

//-----------------------------------------------------------------------------
MeasureDomain::Type MeasureDomain::repr_type(repr_t const& repr)
{
  return MappingToType().find(repr)->second;
}

//-----------------------------------------------------------------------------
MeasureDomain::Type MeasureDomain::type() const
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
                 MeasureData const& meta_data, dolfin::uint& measure_id) :
    Class("Measure"),
    measure_domain_(measure_type),
    measure_id_(measure_id),
    meta_data_(meta_data),
    repr_(*this, measure_domain_, measure_id_, meta_data_),
    str_(measure_domain_.str() + measure_id_.str() + meta_data_.str())
{
}

//-----------------------------------------------------------------------------
Measure::Measure(repr_t const & repr) :
    Class("Measure", repr),
    measure_domain_(arg(0)),
    measure_id_(arg(1)),
    meta_data_(arg(2)),
    repr_(*this, measure_domain_, measure_id_, meta_data_),
    str_(measure_domain_.str() + measure_id_.str() + meta_data_.str())
{
}

//-----------------------------------------------------------------------------
Measure::~Measure()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Measure::operands(std::string const& name) const
{
  std::vector<Class const*> obj0;
  if (name == this->name()) obj0.push_back(this);
  return obj0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Measure::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0 = operands;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);
  return new_operands0;
}

//-----------------------------------------------------------------------------
Measure const* Measure::create(Object::repr_t const& repr)
{
  return new Measure(repr);
}

//-----------------------------------------------------------------------------
MeasureDomain::Type Measure::measure_type() const
{
  return measure_domain_.type();
}

//-----------------------------------------------------------------------------
dolfin::uint Measure::measure_domain_id() const
{
  return measure_id_;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const Measure::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    dolfin::UFCReferenceCell const& ref_cell,
    std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Measure::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Measure::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Measure::display() const
{
}

//-----------------------------------------------------------------------------
Integral::Integral(Expression const& integrand, Measure const& measure) :
    Class("Integral"),
    expressions_(fill_expressions(integrand)),
    measure_(measure),
    repr_(*this, integrand, measure),
    str_("{ " + integrand.str() + " } * " + measure.str())
{
}

//-----------------------------------------------------------------------------
Integral::Integral(repr_t const & repr) :
    Class("Integral", repr),
    expressions_(fill_expressions(args())),
    measure_(arg(1)),
    repr_(*this, *expressions_[0], measure_),
    str_("{ " + expressions_[0]->str() + " } * " + measure_.str())
{
}

//-----------------------------------------------------------------------------
Integral::~Integral()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Integral::operands(
    std::string const& name) const
{
  std::vector<Class const*> obj0 = expressions_[0]->operands(name);
  if (name == this->name()) obj0.push_back(this);
  return obj0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Integral::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0 =
      expressions_[0]->level_operands(operands);
  std::vector<std::vector<Class const*> > new_operands1 =
      measure_.level_operands(operands);
  
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
Integral const* Integral::create(Object::repr_t const& repr)
{
  return new Integral(repr);
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const& Integral::integrand() const
{
  return expressions_;
}

//-----------------------------------------------------------------------------
Measure const& Integral::measure() const
{
  return measure_;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const Integral::evaluate(
    dolfin::uint n,
    std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
    dolfin::UFCReferenceCell const& ref_cell,
    std::vector<dolfin::real*> const& q_points,
    const double * const * coordinates) const
{
  expressions_[0]->evaluate(n, tensor, ref_cell, q_points, coordinates);
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Integral::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Integral::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Integral::display() const
{
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> Integral::fill_expressions(
    std::vector<repr_t> const& reprs)
{
  std::vector<Expression const *> expr;
  expr.push_back(Expression::create(reprs[0]));
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> Integral::fill_expressions(Expression const& e)
{
  std::vector<Expression const *> expr;
  expr.push_back(&e);
  return expr;
}

}
