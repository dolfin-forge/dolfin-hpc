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

#include <dolfin/fem/Equation.h>

#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Equation::Equation() :
    a_(NULL),
    L_(NULL)
{
}

//-----------------------------------------------------------------------------
Equation::~Equation()
{
}

//-----------------------------------------------------------------------------
void Equation::assemble(Matrix& A, Vector& b, bool reset_tensor)
{
  dolfin_assert(a_ != NULL);
  a_->assemble(A, reset_tensor);
  dolfin_assert(L_ != NULL);
  L_->assemble(b, reset_tensor);
}

//-----------------------------------------------------------------------------
bool Equation::is_initialized() const
{
  return (this->a_ != NULL) && (this->L_ != NULL);
}

//-----------------------------------------------------------------------------
void Equation::disp() const
{
  section("Equation");
  section("Bilinear form");
  if (this->a_ != NULL)
  {
    for (uint i = 0; i < a_->num_coefficients(); ++i)
    {
      message("Coefficient %2d : %s", i, this->a_->coefficient_name(i).c_str());
    }
  }
  else
  {
    message("Empty");
  }
  end();
  section("Linear form");
  if (this->L_ != NULL)
  {
    for (uint i = 0; i < this->L_->num_coefficients(); ++i)
    {
      message("Coefficient %2d : %s", i, this->L_->coefficient_name(i).c_str());
    }
  }
  else
  {
    message("Empty");
  }
  end();
  end();
}

//-----------------------------------------------------------------------------
void Equation::clear()
{
  delete a_;
  a_ = NULL;
  delete L_;
  L_ = NULL;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

