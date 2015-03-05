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
  while (-1 != (flag = getopt(argc, argv, "d:m:b")))
  {
    switch (flag)
      {
      case 'b':
        args.benchmark = true;
        break;
      case 'd':
        LogManager::logger().setDebugLevel(std::atoi(optarg));
        break;
      case 'm':
        args.mesh_file = optarg;
        break;
      default:
        message("Unknown or missing argument");
        dolfin_finalize();
        std::exit(1);
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
    std::cout << "mesh : " << args.mesh_file << std::endl;
  }
}

//-----------------------------------------------------------------------------
void Test::begin_test(std::string const& description)
{
  begin("TEST:");
  message(description);
  tic();
}

//-----------------------------------------------------------------------------
void Test::end_test()
{
  message("Completed in %16f s.", toc());
  end();
  skip();
}

//-----------------------------------------------------------------------------
Test::~Test()
{
  dolfin_finalize();
}

}
