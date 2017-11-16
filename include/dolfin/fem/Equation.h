/******************************************************************************
 * Copyright 2013 Aurélien Larcher
 *
 * Licensed under the EUPL, Version 1.1 only (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * http://ec.europa.eu/idabc/eupl5
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 ******************************************************************************/

#ifndef DOLFIN_EQUATION_H_
#define DOLFIN_EQUATION_H_

#include <dolfin/common/types.h>

#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/LinearForm.h>

namespace dolfin
{

class Matrix;
class Vector;

struct Equation
{

  //--- PUBLIC ATTRIBUTES -----------------------------------------------------

  // Bilinear form [Stable]
  BilinearForm * a;

  // Linear form [Stable]
  LinearForm * L;

  //---------------------------------------------------------------------------

  //
  Equation();

  //
  virtual ~Equation();

  //
  void assemble(Matrix& A, Vector& b, bool reset_tensor);

  //
  bool is_initialized() const;

  //
  void disp() const;

  //
  void clear();

  //---------------------------------------------------------------------------

  struct None
  {
    struct BilinearForm : dolfin::BilinearForm
    {
      BilinearForm(Mesh& mesh, CoefficientMap& coefs) :
        dolfin::BilinearForm(mesh) { error("BilinearForm: undefined"); }
      Array<Coefficient*> const& coefficients() const { return c_; }
      ufc::form const& form() const { return *this; }
      Array<Coefficient*> c_;
    };
    struct LinearForm : dolfin::LinearForm
    {
      LinearForm(Mesh& mesh, CoefficientMap& coefs) :
        dolfin::LinearForm(mesh) { error("LinearForm: undefined"); }
      Array<Coefficient*> const& coefficients() const { return c_; }
      ufc::form const& form() const { return *this; }
      Array<Coefficient*> c_;
    };


  };

};

//-----------------------------------------------------------------------------
// Stuck in C++98-land

class CoefficientMap;

template<class E1, class E2, class E3>
struct Equations : public Equation
{

  void operator()(Mesh& mesh, CoefficientMap& coefs)
  {
    Equations<E1, E2, E3>::init(*this, mesh, coefs);
  }

  static inline
  void init(Equation& E, Mesh& mesh, CoefficientMap& coefs)
  {
    if(E.is_initialized())
    {
      error("Equation: already initialized");
    }
    switch (mesh.topology().dim())
    {
      case 1:
        E.a = new typename E1::BilinearForm(mesh, coefs);
        E.L = new typename E1::LinearForm  (mesh, coefs);
        break;
      case 2:
        E.a = new typename E2::BilinearForm(mesh, coefs);
        E.L = new typename E2::LinearForm  (mesh, coefs);
        break;
      case 3:
        E.a = new typename E3::BilinearForm(mesh, coefs);
        E.L = new typename E3::LinearForm  (mesh, coefs);
        break;
      default:
        error("Equation: invalid topological dimension");
        break;
    }
  }
};

} /* namespace dolfin */

#endif /* DOLFIN_EQUATION_H_ */
