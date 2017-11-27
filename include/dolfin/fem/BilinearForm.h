// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#ifndef __DOLFIN_BILINEAR_FORM_H
#define __DOLFIN_BILINEAR_FORM_H

#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/CoefficientMap.h>

#include <ufc.h>

namespace dolfin
{

class GenericMatrix;
class GenericVector;

class BilinearForm : public Form
{
public:

  typedef Form::Coefficients Coefficients;

  static inline std::string name() { return "BilinearForm"; }

  /// Constructor
  BilinearForm(Mesh& mesh);

  /// Destructor
  ~BilinearForm();

  /// Trial space
  FiniteElementSpace const& trial_space() const;

  /// Test space
  FiniteElementSpace const& test_space() const;

  /// Check whether linear system's dimensions match discrete spaces
  void check(GenericMatrix const& A, GenericVector const& b) const;

  /// Creator function
  template <class E> static inline
  typename E::BilinearForm * create(Mesh& mesh, Coefficients& coefs)
  {
    return new typename E::BilinearForm(mesh, coefs);
  }

private:

  mutable FiniteElementSpace * test_space_;
  mutable FiniteElementSpace * trial_space_;

};

//-----------------------------------------------------------------------------

struct NoBilinearForm
{
  typedef Nil<dolfin::BilinearForm> BilinearForm;
};

//--- INLINES -----------------------------------------------------------------

inline FiniteElementSpace const& BilinearForm::trial_space() const
{
  if (!trial_space_)
  {
    trial_space_ = this->create_space(1);
  }
  return *trial_space_;
}

//-----------------------------------------------------------------------------
inline FiniteElementSpace const& BilinearForm::test_space() const
{
  if (!test_space_)
  {
    test_space_ = this->create_space(0);
  }
  return *test_space_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
