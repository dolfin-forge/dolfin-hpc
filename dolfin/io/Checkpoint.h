// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-08
// Last changed: 2009-09-08

#ifndef __CHECKPOINT_H
#define __CHECKPOINT_H

#include <fstream>
#include <string>
#include <vector>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/function/Function.h>

namespace dolfin
{  
  class Checkpoint 
  {
  public:
    
    Checkpoint();
    ~Checkpoint();
    

    void write(real t, Mesh& mesh, std::vector<Function *> func);

    void restart(std::string fname);

    void load(Mesh& mesh);
    void load(std::vector<Function *> func);

  private:

    enum CheckpointState {OPEN, MESH, FUNC};

    CheckpointState state;
    
    std::ifstream in;
  };
}
#endif
