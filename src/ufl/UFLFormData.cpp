// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLFormData.h>
namespace ufl
{

//-----------------------------------------------------------------------------
FormData::FormData(Form const& form) :
    cell_integrals_(form.cell_integrals().operands()),
    extf_integrals_(form.exterior_facet_integrals().operands()),
    intf_integrals_(form.interior_facet_integrals().operands()),
    cell_integrands_(init_cell_integrands()),
    extf_integrands_(init_extf_integrands()),
    intf_integrands_(init_intf_integrands()),
    cell_arguments_(init_cell_arguments()),
    extf_arguments_(init_extf_arguments()),
    intf_arguments_(init_intf_arguments()),
    cell_elements_(init_cell_elements()),
    extf_elements_(init_extf_elements()),
    intf_elements_(init_intf_elements()),
    cell_coefficients_(init_cell_coefficients()),
    extf_coefficients_(init_extf_coefficients()),
    intf_coefficients_(init_intf_coefficients()),
    cell_operands_(init_cell_operands()),
    extf_operands_(init_extf_operands()),
    intf_operands_(init_intf_operands())
{
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::rank() const
{
  std::cout << "cell arg = " << cell_arguments_.size() << std::endl;
  if (cell_arguments_.size() != 0)
  {
    dolfin::uint max_count = 0;
    for (dolfin::uint i = 0; i < cell_arguments_.size(); ++i)
      if (cell_arguments_[i]->count() > max_count) max_count =
          cell_arguments_[i]->count();
    std::cout << "max count = " << max_count << std::endl;
    return max_count + 1;
  }

  std::cout << "extf el = " << extf_arguments_.size() << std::endl;
  if (extf_arguments_.size() != 0)
  {
    dolfin::uint max_count = 0;
    for (dolfin::uint i = 0; i < extf_arguments_.size(); ++i)
      if (extf_arguments_[i]->count() > max_count) max_count =
          extf_arguments_[i]->count();
    std::cout << "max count = " << max_count << std::endl;
    return max_count + 1;
  }

  std::cout << "intf el = " << intf_arguments_.size() << std::endl;
  if (intf_arguments_.size() != 0)
  {
    dolfin::uint max_count = 0;
    for (dolfin::uint i = 0; i < intf_arguments_.size(); ++i)
      if (intf_arguments_[i]->count() > max_count) max_count =
          intf_arguments_[i]->count();
    std::cout << "max count = " << max_count << std::endl;
    return max_count + 1;
  }

  return 0;
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_elements() const
{
  return cell_elements_.size() + extf_elements_.size() + intf_elements_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_cell_elements() const
{
  return cell_elements_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_extf_elements() const
{
  return extf_elements_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_intf_elements() const
{
  return intf_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementSpace const * FormData::get_element(dolfin::uint i) const
{
  std::cout << "i= " << i << " size= " << cell_elements_.size() << std::endl;
  if (cell_elements_.size() > i) if (cell_elements_[i] != NULL) return cell_elements_[i];

  if (extf_elements_.size() > i) if (extf_elements_[i] != NULL) return extf_elements_[i];

  if (intf_elements_.size() > i) if (intf_elements_[i] != NULL) return intf_elements_[i];

  dolfin::error("Either no element present or index %d is out of range", i);
  return NULL;
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::local_finite_element_dimension(dolfin::uint i) const
{
  dolfin::FiniteElement fe(*get_element(i));
  return fe.space_dimension();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_coefficients() const
{
  std::cout << "n cell coefs=" << cell_coefficients_.size() << std::endl;
  std::cout << "n extf coefs=" << extf_coefficients_.size() << std::endl;
  std::cout << "n intf coefs=" << intf_coefficients_.size() << std::endl;
  return cell_coefficients_.size() + extf_coefficients_.size()
      + intf_coefficients_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_cell_coefficients() const
{
  return cell_coefficients_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_extf_coefficients() const
{
  return extf_coefficients_.size();
}

//-----------------------------------------------------------------------------
dolfin::uint FormData::n_intf_coefficients() const
{
  return intf_coefficients_.size();
}

//-----------------------------------------------------------------------------
CoefficientBase const * FormData::get_coefficient(dolfin::uint i) const
{
  std::cout << "i= " << i << " size= " << cell_coefficients_.size()
            << std::endl;
  if (cell_coefficients_.size() > i) if (cell_coefficients_[i] != NULL) return cell_coefficients_[i];

  if (extf_coefficients_.size() > i) if (extf_coefficients_[i] != NULL) return extf_coefficients_[i];

  if (intf_coefficients_.size() > i) if (intf_coefficients_[i] != NULL) return intf_coefficients_[i];

  dolfin::error("Either no element present or index %d is out of range", i);
  return NULL;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const *> > const& FormData::get_cell_operands(
    dolfin::uint i) const
{
  return cell_operands_[i];
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const *> > const& FormData::get_extf_operands(
    dolfin::uint i) const
{
  return extf_operands_[i];
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const *> > const& FormData::get_intf_operands(
    dolfin::uint i) const
{
  return intf_operands_[i];
}

//-----------------------------------------------------------------------------
void FormData::zero_element_tensor(double* A) const
{
  unsigned int tensor_rank = this->rank();
  std::vector<unsigned int> n_test(tensor_rank);
  for (unsigned int i = 0; i < n_test.size(); ++i)
    n_test[i] = local_finite_element_dimension(i);

  for (unsigned int j = 0; j < n_test[0]; j++)
  {
    if (tensor_rank == 2)
    {
      for (unsigned int k = 0; k < n_test[1]; k++)
      {
        A[j * n_test[0] + k] = 0;
      }  // end loop over 'k'
    }
    else if (tensor_rank == 1) A[j] = 0;
    else std::cout << "not implemented" << std::endl;

  }  // end loop over 'j'
}

//-----------------------------------------------------------------------------
void FormData::compute_element_tensor(double*,
                                      dolfin::UFCReferenceCell const& ref_cell,
                                      std::vector<dolfin::real*> const& q_point,
                                      std::vector<dolfin::real> const&,
                                      const double * const * coordinates) const
{
  std::vector<std::vector<std::vector<std::vector<dolfin::real> > > > new_vals0(
      cell_integrals_.size());
  for (dolfin::uint i = 0; i < cell_integrals_.size(); ++i)
  {
    new_vals0[i] = cell_integrals_[i]->evaluate(0, new_vals0[i], ref_cell,
                                                q_point, coordinates);

//        for(dolfin::uint k=0; k<expr[i].size(); ++k)
//          for(dolfin::uint j=0; j<expr[i][k].size(); ++j)
//            std::cout << "(" << k << "," << j << ")" << expr[i][k][j]->name() << std::endl;
  }
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::get_cell_arguments() const
{
  return cell_arguments_;
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::get_extf_arguments() const
{
  return extf_arguments_;
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::get_intf_arguments() const
{
  return intf_arguments_;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::get_cell_elements() const
{
  return cell_elements_;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::get_extf_elements() const
{
  return extf_elements_;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::get_intf_elements() const
{
  return intf_elements_;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::get_cell_coefficients() const
{
  return cell_coefficients_;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::get_extf_coefficients() const
{
  return extf_coefficients_;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::get_intf_coefficients() const
{
  return intf_coefficients_;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const FormData::init_cell_integrands() const
{
  std::vector<Expression const *> cell_integrands;
  for (dolfin::uint i = 0; i < cell_integrals_.size(); ++i)
  {
    std::vector<Expression const *> const& integrands =
        cell_integrals_[i]->integrand();
    for (dolfin::uint j = 0; j < integrands.size(); ++j)
      cell_integrands.push_back(integrands[j]);
  }

  return cell_integrands;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const FormData::init_extf_integrands() const
{
  std::vector<Expression const *> extf_integrands;
  for (dolfin::uint i = 0; i < extf_integrals_.size(); ++i)
  {
    std::vector<Expression const *> const& integrands =
        extf_integrals_[i]->integrand();
    for (dolfin::uint j = 0; j < integrands.size(); ++j)
      extf_integrands.push_back(integrands[j]);
  }

  return extf_integrands;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const FormData::init_intf_integrands() const
{
  std::vector<Expression const *> intf_integrands;
  for (dolfin::uint i = 0; i < intf_integrals_.size(); ++i)
  {
    std::vector<Expression const *> const& integrands =
        intf_integrals_[i]->integrand();
    for (dolfin::uint j = 0; j < integrands.size(); ++j)
      intf_integrands.push_back(integrands[j]);
  }

  return intf_integrands;
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::init_cell_arguments() const
{
  std::vector<Argument const *> cell_arguments;
  for (dolfin::uint i = 0; i < cell_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& arguments = cell_integrands_[i]->operands(
        "Argument");
    for (dolfin::uint j = 0; j < arguments.size(); ++j)
      cell_arguments.push_back(dynamic_cast<Argument const *>(arguments[j]));
  }

  return cell_arguments;
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::init_extf_arguments() const
{
  std::vector<Argument const *> extf_arguments;
  for (dolfin::uint i = 0; i < extf_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& arguments = extf_integrands_[i]->operands(
        "Argument");
    for (dolfin::uint j = 0; j < arguments.size(); ++j)
      extf_arguments.push_back(dynamic_cast<Argument const *>(arguments[j]));
  }

  return extf_arguments;
}

//-----------------------------------------------------------------------------
std::vector<Argument const *> const FormData::init_intf_arguments() const
{
  std::vector<Argument const *> intf_arguments;
  for (dolfin::uint i = 0; i < intf_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& arguments = intf_integrands_[i]->operands(
        "Argument");
    for (dolfin::uint j = 0; j < arguments.size(); ++j)
      intf_arguments.push_back(dynamic_cast<Argument const *>(arguments[j]));
  }

  return intf_arguments;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::init_cell_elements() const
{
  std::vector<FiniteElementSpace const *> cell_elements;
  for (dolfin::uint i = 0; i < cell_arguments_.size(); ++i)
    cell_elements.push_back(&cell_arguments_[i]->element());

  std::cout << "size cell_elements= " << cell_elements.size() << std::endl;
  return cell_elements;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::init_extf_elements() const
{
  std::vector<FiniteElementSpace const *> extf_elements;
  for (dolfin::uint i = 0; i < extf_arguments_.size(); ++i)
    extf_elements.push_back(&extf_arguments_[i]->element());

  return extf_elements;
}

//-----------------------------------------------------------------------------
std::vector<FiniteElementSpace const *> const FormData::init_intf_elements() const
{
  std::vector<FiniteElementSpace const *> intf_elements;
  for (dolfin::uint i = 0; i < intf_arguments_.size(); ++i)
    intf_elements.push_back(&intf_arguments_[i]->element());

  return intf_elements;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::init_cell_coefficients() const
{
  std::vector<CoefficientBase const *> cell_coefficients;
  for (dolfin::uint i = 0; i < cell_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& coef = cell_integrands_[i]->operands(
        "Coefficient");
    for (dolfin::uint j = 0; j < coef.size(); ++j)
      cell_coefficients.push_back(dynamic_cast<Coefficient const *>(coef[j]));
    std::vector<Class const *> const& cnst = cell_integrands_[i]->operands(
        "Constant");
    for (dolfin::uint j = 0; j < cnst.size(); ++j)
      cell_coefficients.push_back(dynamic_cast<Constant const *>(cnst[j]));
    std::vector<Class const *> const& vcnst = cell_integrands_[i]->operands(
        "VectorConstant");
    for (dolfin::uint j = 0; j < vcnst.size(); ++j)
      cell_coefficients.push_back(
          dynamic_cast<VectorConstant const *>(vcnst[j]));
    std::vector<Class const *> const& tcnst = cell_integrands_[i]->operands(
        "TensorConstant");
    for (dolfin::uint j = 0; j < tcnst.size(); ++j)
      cell_coefficients.push_back(
          dynamic_cast<TensorConstant const *>(tcnst[j]));
  }

  return cell_coefficients;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::init_extf_coefficients() const
{
  std::vector<CoefficientBase const *> extf_coefficients;
  for (dolfin::uint i = 0; i < extf_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& coef = extf_integrands_[i]->operands(
        "Coefficient");
    for (dolfin::uint j = 0; j < coef.size(); ++j)
      extf_coefficients.push_back(dynamic_cast<Coefficient const *>(coef[j]));
    std::vector<Class const *> const& cnst = extf_integrands_[i]->operands(
        "Constant");
    for (dolfin::uint j = 0; j < cnst.size(); ++j)
      extf_coefficients.push_back(dynamic_cast<Constant const *>(cnst[j]));
    std::vector<Class const *> const& vcnst = extf_integrands_[i]->operands(
        "VectorConstant");
    for (dolfin::uint j = 0; j < vcnst.size(); ++j)
      extf_coefficients.push_back(
          dynamic_cast<VectorConstant const *>(vcnst[j]));
    std::vector<Class const *> const& tcnst = extf_integrands_[i]->operands(
        "TensorConstant");
    for (dolfin::uint j = 0; j < tcnst.size(); ++j)
      extf_coefficients.push_back(
          dynamic_cast<TensorConstant const *>(tcnst[j]));
  }

  return extf_coefficients;
}

//-----------------------------------------------------------------------------
std::vector<CoefficientBase const *> const FormData::init_intf_coefficients() const
{
  std::vector<CoefficientBase const *> intf_coefficients;
  for (dolfin::uint i = 0; i < intf_integrands_.size(); ++i)
  {
    std::vector<Class const *> const& coef = intf_integrands_[i]->operands(
        "Coefficient");
    for (dolfin::uint j = 0; j < coef.size(); ++j)
      intf_coefficients.push_back(dynamic_cast<Coefficient const *>(coef[j]));
    std::vector<Class const *> const& cnst = intf_integrands_[i]->operands(
        "Constant");
    for (dolfin::uint j = 0; j < cnst.size(); ++j)
      intf_coefficients.push_back(dynamic_cast<Constant const *>(cnst[j]));
    std::vector<Class const *> const& vcnst = intf_integrands_[i]->operands(
        "VectorConstant");
    for (dolfin::uint j = 0; j < vcnst.size(); ++j)
      intf_coefficients.push_back(
          dynamic_cast<VectorConstant const *>(vcnst[j]));
    std::vector<Class const *> const& tcnst = intf_integrands_[i]->operands(
        "TensorConstant");
    for (dolfin::uint j = 0; j < tcnst.size(); ++j)
      intf_coefficients.push_back(
          dynamic_cast<TensorConstant const *>(tcnst[j]));
  }

  return intf_coefficients;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<Class const *> > > const FormData::init_cell_operands() const
{
  std::vector<std::vector<std::vector<Class const *> > > expr(
      cell_integrals_.size());
  for (dolfin::uint i = 0; i < cell_integrals_.size(); ++i)
  {
    expr[i] = cell_integrals_[i]->level_operands(expr[i]);

//        for(dolfin::uint k=0; k<expr[i].size(); ++k)
//          for(dolfin::uint j=0; j<expr[i][k].size(); ++j)
//            std::cout << "(" << k << "," << j << ")" << expr[i][k][j]->name() << std::endl;
  }
  return expr;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<Class const *> > > const FormData::init_extf_operands() const
{
  std::vector<std::vector<std::vector<Class const *> > > expr(
      extf_integrals_.size());
  for (dolfin::uint i = 0; i < extf_integrals_.size(); ++i)
  {
    expr[i] = extf_integrals_[i]->level_operands(expr[i]);

//        for(dolfin::uint k=0; k<expr[i].size(); ++k)
//          for(dolfin::uint j=0; j<expr[i][k].size(); ++j)
//            std::cout << "(" << k << "," << j << ")" << expr[i][k][j]->name() << std::endl;
  }

  return expr;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<Class const *> > > const FormData::init_intf_operands() const
{
  std::vector<std::vector<std::vector<Class const *> > > expr(
      intf_integrals_.size());
  for (dolfin::uint i = 0; i < intf_integrals_.size(); ++i)
  {
    expr[i] = intf_integrals_[i]->level_operands(expr[i]);

//        for(dolfin::uint k=0; k<expr[i].size(); ++k)
//          for(dolfin::uint j=0; j<expr[i][k].size(); ++j)
//            std::cout << "(" << k << "," << j << ")" << expr[i][k][j]->name() << std::endl;
  }

  return expr;
}

}
