// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#ifndef __DOLFIN_TEST_H
#define __DOLFIN_TEST_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

#include <iostream>
#include <string>

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
  void begin(std::string const& name);

  ///
  void end();

  ///
  ~Test();

  Args args;

private:

  bool btest_;
  Array<std::pair<std::string, real> > timings_;
  real total_;
  uint padding_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_TEST_H */
