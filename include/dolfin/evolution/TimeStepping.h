/******************************************************************************
 * Copyright 2012, 2013 Aurélien Larcher
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

#ifndef DOLFIN_TIME_STEPPING_H_
#define DOLFIN_TIME_STEPPING_H_

#include <licorne/common/types.h>

namespace dolfin
{

/**
 * 	DOCUMENTATION:
 *
 * 	@class    TimeStepping
 *
 * 	@brief    Provides an interface to instantiate time stepping schemes.
 *
 *	@author   Aurélien Larcher <larcher@kth.se>
 *	@version  1
 */

class TimeStepping
{

public:

  ///
  TimeStepping() {}

  ///
  virtual ~TimeStepping() {}

  //--- INTERFACE -------------------------------------------------------------

  ///
  virtual void update() = 0;

  ///
  virtual void evolve() = 0;

  ///
  virtual void derive() = 0;

  ///
  virtual uint order() const = 0;

  //---------------------------------------------------------------------------

};

} /* namespace dolfin */

#endif /* DOLFIN_TIME_STEPPING_H_ */
