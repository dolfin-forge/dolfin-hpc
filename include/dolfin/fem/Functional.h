// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTIONAL_H
#define __DOLFIN_FUNCTIONAL_H

#include <dolfin/fem/Form.h>
#include <dolfin/fem/CoefficientMap.h>

#include <ufc.h>

namespace dolfin
{

class Functional : public Form
{

public:

  typedef Form::Coefficients Coefficients;

  static inline std::string name() { return "Functional"; }

  /// Constructor
  Functional(Mesh& mesh);

  /// Destructor
  ~Functional();

  /// Creator function
  template <class E> static inline
  typename E::Functional * create(Mesh& mesh, Coefficients& coefs)
  {
    return new typename E::Functional(mesh, coefs);
  }

private:

};

//-----------------------------------------------------------------------------

struct NoFunctional
{
  typedef Nil<dolfin::Functional> Functional;
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
