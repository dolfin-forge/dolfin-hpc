// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_FUNCTION_INTERPOLATION_H
#define __DOLFIN_FUNCTION_INTERPOLATION_H

namespace dolfin
{

class Coefficient;
class Expression;
class Function;
class GenericFunction;

class FunctionInterpolation
{

public:

  ///
  static void compute(Expression const& F0, Function& F1);

  ///
  static void compute(Coefficient const& F0, Function& F1);

  ///
  static void compute(GenericFunction const& F0, Function& F1);

  ///
  ~FunctionInterpolation();

private:

  ///
  static void interpolateSM(GenericFunction const& F0, Function& F1);

  ///
  static void interpolateNM(GenericFunction const& F0, Function& F1);

};

} /* namespace dolfin */

#endif /* __DOLFIN_FUNCTION_INTERPOLATION_H */
