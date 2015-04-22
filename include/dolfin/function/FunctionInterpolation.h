// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __FUNCTION_INTERPOLATION_H_
#define __FUNCTION_INTERPOLATION_H_

namespace dolfin
{

class Function;

class FunctionInterpolation
{

public:

  ///
  FunctionInterpolation(Function const& F0, Function& F1);

  ///
  ~FunctionInterpolation();

  ///
  void compute();

private:

  ///
  void interpolateSameMesh(Function const& F0, Function& F1);

  ///
  void interpolateNonMatchingMeshes(Function const& F0, Function& F1);

  //--- ATTRIBUTES ------------------------------------------------------------

  Function const& F0_;
  Function& F1_;

};

} /* namespace dolfin */

#endif /* __FUNCTION_INTERPOLATION_H_ */
