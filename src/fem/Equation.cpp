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

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Equation::Equation() :
    a(NULL),
    L(NULL)
{
}

//-----------------------------------------------------------------------------
Equation::~Equation()
{
}

//-----------------------------------------------------------------------------
bool Equation::is_initialized() const
{
  return (this->a != NULL) && (this->L != NULL);
}

//-----------------------------------------------------------------------------
void Equation::disp() const
{
  section("Equation");
  section("Bilinear form");
  if (this->a != NULL)
  {
    for (uint i = 0; i < a->coefficients().size(); ++i)
    {
      message("Coefficient %2d : %s", i, this->a->coefficient_name(i).c_str());
    }
  }
  else
  {
    message("Empty");
  }
  end();
  section("Linear form");
  if (this->L != NULL)
  {
    for (uint i = 0; i < this->L->coefficients().size(); ++i)
    {
      message("Coefficient %2d : %s", i, this->L->coefficient_name(i).c_str());
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
  delete a;
  a = NULL;
  delete L;
  L = NULL;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

