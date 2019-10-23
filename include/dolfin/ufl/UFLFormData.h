// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_FORM_DATA_H
#define __DOLFIN_UFL_FORM_DATA_H

#include <dolfin/ufl/UFLArgument.h>
#include <dolfin/ufl/UFLCoefficient.h>
#include <dolfin/ufl/UFLForm.h>
#include <dolfin/quadrature/UFCReferenceCell.h>

namespace ufl
{

class Form;
class UFCReferenceCell;

/**
 *  DOCUMENTATION:
 *
 *  @class  FormData
 *
 *  @brief  Access to data of a ufl Form.
 */

class FormData  // : public Class
{
public:

  FormData(Form const& form);

  dolfin::uint rank() const;

  dolfin::uint n_elements() const;
  dolfin::uint n_cell_elements() const;
  dolfin::uint n_extf_elements() const;
  dolfin::uint n_intf_elements() const;

  FiniteElementSpace const * get_element(dolfin::uint i) const;

  dolfin::uint local_finite_element_dimension(dolfin::uint i) const;

  dolfin::uint n_coefficients() const;
  dolfin::uint n_cell_coefficients() const;
  dolfin::uint n_extf_coefficients() const;
  dolfin::uint n_intf_coefficients() const;

  CoefficientBase const * get_coefficient(dolfin::uint i) const;

  std::vector<std::vector<Class const *> > const& get_cell_operands(
      dolfin::uint i) const;
  std::vector<std::vector<Class const *> > const& get_extf_operands(
      dolfin::uint i) const;
  std::vector<std::vector<Class const *> > const& get_intf_operands(
      dolfin::uint i) const;

  void zero_element_tensor(double* A) const;

  void compute_element_tensor(double* A,
                              dolfin::UFCReferenceCell const& ref_cell,
                              std::vector<dolfin::real*> const& q_point,
                              std::vector<dolfin::real> const& weights,
                              const double * const * coordinates) const;

protected:

  std::vector<Argument const *> const get_cell_arguments() const;
  std::vector<Argument const *> const get_extf_arguments() const;
  std::vector<Argument const *> const get_intf_arguments() const;

  std::vector<FiniteElementSpace const *> const get_cell_elements() const;
  std::vector<FiniteElementSpace const *> const get_extf_elements() const;
  std::vector<FiniteElementSpace const *> const get_intf_elements() const;

  std::vector<CoefficientBase const *> const get_cell_coefficients() const;
  std::vector<CoefficientBase const *> const get_extf_coefficients() const;
  std::vector<CoefficientBase const *> const get_intf_coefficients() const;

private:

  Form const& form_;
  std::vector<Integral const *> const cell_integrals_;
  std::vector<Integral const *> const extf_integrals_;
  std::vector<Integral const *> const intf_integrals_;

  std::vector<Expression const *> const cell_integrands_;
  std::vector<Expression const *> const extf_integrands_;
  std::vector<Expression const *> const intf_integrands_;

  std::vector<Argument const *> const cell_arguments_;
  std::vector<Argument const *> const extf_arguments_;
  std::vector<Argument const *> const intf_arguments_;

  std::vector<FiniteElementSpace const *> const cell_elements_;
  std::vector<FiniteElementSpace const *> const extf_elements_;
  std::vector<FiniteElementSpace const *> const intf_elements_;

  std::vector<CoefficientBase const *> const cell_coefficients_;
  std::vector<CoefficientBase const *> const extf_coefficients_;
  std::vector<CoefficientBase const *> const intf_coefficients_;

  std::vector<std::vector<std::vector<Class const *> > > const cell_operands_;
  std::vector<std::vector<std::vector<Class const *> > > const extf_operands_;
  std::vector<std::vector<std::vector<Class const *> > > const intf_operands_;

  std::vector<Expression const *> const init_cell_integrands() const;
  std::vector<Expression const *> const init_extf_integrands() const;
  std::vector<Expression const *> const init_intf_integrands() const;

  std::vector<Argument const *> const init_cell_arguments() const;
  std::vector<Argument const *> const init_extf_arguments() const;
  std::vector<Argument const *> const init_intf_arguments() const;

  std::vector<FiniteElementSpace const *> const init_cell_elements() const;
  std::vector<FiniteElementSpace const *> const init_extf_elements() const;
  std::vector<FiniteElementSpace const *> const init_intf_elements() const;

  std::vector<CoefficientBase const *> const init_cell_coefficients() const;
  std::vector<CoefficientBase const *> const init_extf_coefficients() const;
  std::vector<CoefficientBase const *> const init_intf_coefficients() const;

  std::vector<std::vector<std::vector<Class const *> > > const init_cell_operands() const;
  std::vector<std::vector<std::vector<Class const *> > > const init_extf_operands() const;
  std::vector<std::vector<std::vector<Class const *> > > const init_intf_operands() const;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_FORM_DATA_H */
