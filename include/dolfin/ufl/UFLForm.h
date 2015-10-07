// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_UFL_FORM_H
#define __DOLFIN_UFL_FORM_H

#include <dolfin/ufl/UFLFormData.h>
#include <dolfin/ufl/UFLList.h>
#include <dolfin/ufl/UFLtuple.h>

#include <ufc.h>

namespace ufl
{

class FormData;

/**
 *  DOCUMENTATION:
 *
 *  @class  Form
 *
 *  @brief  Provides an interface complying with UFL Form.
 */

class Form : public Class
{
public:

  ///
  Form(List const& list);

  ///
  Form(repr_t const& repr);

  ///
  Form(ufc::form const& ufc_form);

  ///
  ~Form();

  //--- INTERFACE -------------------------------------------------------------

  ///Does this function make sense?
  Cell cell() const;

  tuple<Integral> const integrals(MeasureDomain::Type const& measure_type =
      MeasureDomain::None) const;

  tuple<Measure> const measures(MeasureDomain::Type const& measure_type =
      MeasureDomain::None) const;

  std::vector<std::pair<MeasureDomain::Type, dolfin::uint> > const domains(
      MeasureDomain::Type const& measure_type = MeasureDomain::None) const;

  tuple<Integral> const cell_integrals() const;

  tuple<Integral> const exterior_facet_integrals() const;

  tuple<Integral> const interior_facet_integrals() const;

  tuple<Integral> const macro_cell_integrals() const;

  tuple<Integral> const surface_integrals() const;

  FormData const& form_data() const;

//      FormData const& compute_form_data() const;

//      bool const is_preprocessed() const;

//--- INTERFACE inherited from UFLClass -------------------------------------

  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  List const list_;

  FormData const * form_data_;

  repr_t const repr_;
  std::string const str_;

//      bool const is_preprocessed_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_FORM_H */
