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

#ifndef __LICORNE_FIELD_H_
#define __LICORNE_FIELD_H_

#include <dolfin/common/types.h>
#include <dolfin/fem/SetOfBCs.h>

#include <dolfin/common/Array.h>
#include <dolfin/function/Function.h>

namespace dolfin
{

class FiniteElementSpace;

/**
 *
 */

class Field
{

public:

  ///
  Field(Mesh& mesh);

  ///
  Field(FiniteElementSpace const& space);

  ///
  ~Field();

  ///
  Function& dt();

  ///
  uint depth() const;

  ///
  void allocate(uint max_depth);

  ///
  void init(FiniteElementSpace const& space);

  /// Circular shift
  void shift();

  /// Circular n-shift
  Field& operator<<(uint n);

  ///
  Function& operator[](uint i);

  ///
  Function const& operator[](uint i) const;

  ///
  SetOfBCs& bcs();

  ///
  SetOfBCs const& bcs() const;

  ///
  void disp() const;

private:

  Array<Function> storage_;
  SetOfBCs bcs_;
  Function Udt_;

};

} /* namespace dolfin */

#endif /* __LICORNE_FIELD_H_ */
