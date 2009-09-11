// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-08
// Last changed: 2009-09-08

#ifndef __CHECKPOINT_H
#define __CHECKPOINT_H

#include <string>
#include <vector>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/function/Function.h>

namespace dolfin
{  
  class Checkpoint 
  {
  public:

    static void write(real t, Mesh& mesh, std::vector<Function *> func);

    static void restart(std::string fname, 
			Mesh& mesh, std::vector<Function *> func);
  };
}
#endif
