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
  typedef CoefficientMap Coefficients;

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

  //
  void swap(Equation& other)
  {
    std::swap(a_, other.a_);
    std::swap(L_, other.L_);
  }

  //
  inline BilinearForm& a()
  {
    dolfin_assert(a_);
    return *a_;
  }

  //
  inline LinearForm& L()
  {
    dolfin_assert(L_);
    return *L_;
  }

  //---------------------------------------------------------------------------

  struct None
  {
    typedef Nil<dolfin::BilinearForm> BilinearForm;
    typedef Nil<dolfin::LinearForm>   LinearForm;
  };

  /// Creator function
  template <class E>
  void operator()(Mesh& mesh, Coefficients& coefs)
  {
    if (this->is_initialized())
    {
      error("Equation : initialization on non-empty instance");
    }
    a_ = BilinearForm::create<E>(mesh, coefs);
    L_ = LinearForm::create<E>  (mesh, coefs);
  }

  /// Creator function
  template<class E1, class E2, class E3>
  void operator()(Mesh& mesh, Coefficients& coefs)
  {
    switch (mesh.topology().dim())
    {
      case 1:
        Equation::template operator()<E1>(mesh, coefs);
        break;
      case 2:
        Equation::template operator()<E2>(mesh, coefs);
        break;
      case 3:
        Equation::template operator()<E3>(mesh, coefs);
        break;
      default:
        error("Equation : invalid topological dimension");
        break;
    }
  }

  /// Creator function
  template<class E>
  void instantiate(Mesh& mesh, Coefficients& coefs)
  {
    E X; X(mesh, coefs); this->swap(X);
  }

private:

  Equation(Equation const& other) {}

  // Bilinear form [Stable]
  BilinearForm * a_;

  // Linear form [Stable]
  LinearForm * L_;

};

//-----------------------------------------------------------------------------
// Stuck in C++98-land

template<class E1, class E2, class E3>
struct Equations : public Equation
{
  void operator()(Mesh& mesh, Coefficients& coefs)
  {
    Equation::template operator()<E1, E2, E3>(mesh, coefs);
  }
};

} /* namespace dolfin */

#endif /* DOLFIN_EQUATION_H_ */
