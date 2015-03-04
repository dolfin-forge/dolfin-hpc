// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#ifndef __TEST_H
#define __TEST_H

#include <dolfin.h>

#include <cassert>
#include <iostream>

namespace dolfin
{

class Test
{

  struct Args
  {
    uint debug_level;
    std::string mesh_file;
    bool benchmark;

    Args() :
        debug_level(0),
        mesh_file(""),
        benchmark(false)
    {
    }
  };

public:

  ///
  Test(int argc, char *argv[]);

  ///
  Test();

  ///
  void print_args();

  ///
  void begin_test(std::string const& description);

  ///
  void end_test();

  ///
  ~Test();

  Args args;

};

}

#endif /* __TEST_H */
