// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#include <dolfin/common/Test.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogManager.h>
#include <dolfin/main/init.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

#include <cstdlib>
#include <getopt.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Test::Test(int argc, char *argv[])
{
  dolfin_init(argc, argv);
  if (dolfin::MPI::processNumber() == 0)
  {
    message("Running on %d %s", dolfin::MPI::numProcesses(),
            (dolfin::MPI::numProcesses() > 1 ? "nodes" : "node"));
  }

  int flag;
  int i = 0;
  while (-1 != (flag = getopt(argc, argv, ":d:m:b")))
  {
    switch (flag)
      {
      case 'b':
        args.benchmark = true;
        break;
      case 'd':
        args.debug_level = std::atoi(optarg);
        LogManager::logger().setDebugLevel(args.debug_level);
        break;
      case 'm':
        args.mesh_file = optarg;
        break;
      default:
        if(i < argc && (getopt(argc, argv, ":d:m:b") == -1))
        {
          ++optind;
        }
        break;
      }
  }
}

//-----------------------------------------------------------------------------
Test::Test()
{
}

//-----------------------------------------------------------------------------
void Test::print_args()
{
  if (dolfin::MPI::processNumber() == 0)
  {
    message("benchmark   : %s", args.benchmark);
    message("debug level : %u", args.debug_level);
    message("mesh        : %u", args.mesh_file.c_str());
  }
}

//-----------------------------------------------------------------------------
void Test::begin(std::string const& name)
{
  dolfin::begin("TEST:");
  header(name);
  tic();
}

//-----------------------------------------------------------------------------
void Test::end()
{
  tocd();
  dolfin::end();
  skip();
}

//-----------------------------------------------------------------------------
Test::~Test()
{
  dolfin_finalize();
}

}
