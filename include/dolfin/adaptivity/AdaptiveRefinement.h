// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-09-13
// Last changed: 2011-01-18

#ifndef __ADAPTIVEREFINEMENT_H
#define __ADAPTIVEREFINEMENT_H

#include <dolfin/common/types.h>
#include <vector>

namespace dolfin
{

template<class T> class Array;
class Form;
class Function;
class Mesh;
class MeshData;
template<class T> class MeshFunction;
class Function;
class Vector;

class AdaptiveRefinement
{
public:

  ///
  typedef std::pair<Form *, uint> form_tuple;
  typedef std::pair<Function *, form_tuple> project_func;

  ///
  static void refine(Mesh& mesh, MeshFunction<bool>& cell_marker);

  ///
  static void refine_and_project(Mesh& mesh, Array<Function *> const& functions,
                                 MeshFunction<bool>& cell_marker);

private:

  ///
  static void redistribute_func(Mesh& mesh, Function const& f, real **vp,
                                uint **rp, uint& m,
                                MeshFunction<uint>& distribution);

  ///
  static void decompose_func(Mesh& mesh, Function const& function,
                             Array<Function *>& subfunctions);

  ///
  static void project(Mesh& new_mesh, Array<Function *>& f_post,
                      Function& projected);

};
}
#endif
