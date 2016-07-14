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
  bool is_initialized() const;

  //
  void disp() const;

  //
  void clear();

};

} /* namespace dolfin */

#endif /* DOLFIN_EQUATION_H_ */
