// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TEST_H
#define __DOLFIN_TEST_H

#include <dolfin/main/Startup.h>

#include <dolfin/common/Array.h>

#include <iostream>
#include <string>

namespace dolfin
{

class Test : public Startup
{

  struct Args
  {
    uint debug_level{0};
    std::string mesh_file;
    bool benchmark{false};

    Args()
      : mesh_file("")
    {
    }
  };

public:

  ///
  Test(int argc, char *argv[]);

  ///
  Test(std::string const& dir = "");

  ///
  void print_args();

  ///
  void begin(std::string const& name);

  ///
  void end();

  ///
  ~Test() override;

  Args args;

private:

  ///
  void init(int argc, char *argv[]);

  //--- ATTRIBUTES ------------------------------------------------------------

  bool btest_{false};
  std::string const dir_;
  Array<std::pair<std::string, real> > timings_;
  real total_{0.0};
  uint padding_{0};

};

} /* namespace dolfin */

#endif /* __DOLFIN_TEST_H */
