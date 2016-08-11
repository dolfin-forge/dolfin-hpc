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

#ifndef LICORNE_SET_OF_BCS_H_
#define LICORNE_SET_OF_BCS_H_

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/BoundaryCondition.h>

namespace dolfin
{

class BilinearForm;
class GenericMatrix;
class GenericVector;

/**
 *  DOCUMENTATION:
 *
 *  @class SetOfBCs
 *
 *  @brief
 *
 *  @author   Aurélien Larcher <larcher@kth.se>
 *  @version  1
 */

class SetOfBCs : private Array<BoundaryCondition *>
{
  typedef Array<BoundaryCondition *> ArrayOfBCs;

public:

  /// Constructor
  SetOfBCs();

  /// Destructor
  ~SetOfBCs();

  /// Add boundary condition
  void add(BoundaryCondition& bc);

  /// Apply all boundary conditions to the given system
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm& a) const;

  /// Apply all boundary conditions to the given subsystem
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm& a,
             SubSystem const sub) const;

  /// Return if the given boundary condition type exists
  bool has(std::string const& type) const;

  /// Return all the boundary conditions for given type
  Array<BoundaryCondition *> get(std::string const& type) const;

  ///
  uint size() const;

  /// Display the list of boundary conditions
  void disp() const;

private:

  ///
  void clear();

};

} /* namespace dolfin */

#endif /* LICORNE_SET_OF_BCS_H_ */
