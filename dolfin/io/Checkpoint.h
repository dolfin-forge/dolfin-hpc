// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-08
// Last changed: 2009-09-14

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
    

    void write(uint id, real t, Mesh& mesh, std::vector<Function *> func);    
    void restart(std::string fname);

    void load(Mesh& mesh);
    void load(std::vector<Function *> func);

    inline bool restart() {return state == RESTART;};

    inline dolfin::uint id()
    { if(state != RESTART)
	error("Shut her down, Scotty, she's sucking mud again!");
	return _id;
    };

    inline dolfin::real restart_time()
    { if(state != RESTART) 
	error("Shut her down, Scotty, she's sucking mud again!");
      return _t;
    };

    inline void reset() 
    { state = CHECKPOINT; restart_state = OPEN; };

  private:
    
    void write(Mesh& mesh, std::ofstream& out);
    void write(std::vector<Function *> func, std::ofstream& out);

    enum CheckpointState {CHECKPOINT, RESTART};
    enum RestartState {OPEN, MESH, FUNC};

    CheckpointState state;
    RestartState restart_state;
    
    std::ifstream in;
    
    uint _id;
    real _t;
  };
}
#endif
