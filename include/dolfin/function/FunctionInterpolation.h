// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_FUNCTION_INTERPOLATION_H
#define __DOLFIN_FUNCTION_INTERPOLATION_H

namespace dolfin
{

class Function;
class GenericFunction;

class FunctionInterpolation
{

public:

  ///
  explicit FunctionInterpolation(GenericFunction const& F0, Function& F1);

  ///
  ~FunctionInterpolation();

  ///
  void compute();

private:

  ///
  void interpolateSM(GenericFunction const& F0, Function& F1);

  ///
  void interpolateNM(GenericFunction const& F0, Function& F1);

  //--- ATTRIBUTES ------------------------------------------------------------

  GenericFunction const& F0_;
  Function& F1_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_FUNCTION_INTERPOLATION_H */
