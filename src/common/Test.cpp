// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/Test.h>

#include <dolfin/common/timing.h>
#include <dolfin/common/system.h>
#include <dolfin/log/log.h>
#include <dolfin/main/init.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

#include <sstream>
#include <cstdlib>
#include <getopt.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Test::Test(int argc, char *argv[]) :
    Startup(argc, argv),

    dir_("")

{
  init(argc, argv);
}

//-----------------------------------------------------------------------------
Test::Test(std::string const& dir) :
    Startup(0, nullptr),
    dir_(dir)
{
  init(0, nullptr);
  if(!dir_.empty())
  {
    pushd(dir_);
  }
}

//-----------------------------------------------------------------------------
void Test::init(int argc, char *argv[])
{
  int flag;
  while (argc > 1 && -1 != (flag = getopt(argc, argv, "d:m:b")))
  {
    switch (flag)
      {
      case 'b':
        args.benchmark = true;
        break;
      case 'd':
        args.debug_level = std::atoi(optarg);
        verbose(args.debug_level);
        break;
      case 'm':
        args.mesh_file = optarg;
        break;
      default:
        break;
      }
  }
  // Reset optind if applications are using getopt through the TEST interface
  optind = 1;
}

//-----------------------------------------------------------------------------
void Test::print_args()
{
  if (dolfin::MPI::rank() == 0)
  {
    message("benchmark   : %s", args.benchmark);
    message("debug level : %u", args.debug_level);
    message("mesh        : %u", args.mesh_file.c_str());
  }
}

//-----------------------------------------------------------------------------
void Test::begin(std::string const& name)
{
  if(btest_)
  {
    error("Test : missing end directive for preceding test");
  }
  btest_ = true;
  timings_.push_back(std::pair<std::string, real >(name, 0.0));
  padding_ = std::max(padding_, name.size());
  std::stringstream ss;
  ss << "Test " << timings_.size() << " : " << name;
  message(ss.str());
  skip();
  tic();
}

//-----------------------------------------------------------------------------
void Test::end()
{
  if(!btest_)
  {
    error("Test : missing end directive for preceding test");
  }
  btest_ = false;
  real elapsed_time = toc();
  cout << "Elapsed time: " << elapsed_time << " seconds\n";
  timings_.back().second = elapsed_time;
  total_ += elapsed_time;
  skip();
}

//-----------------------------------------------------------------------------
Test::~Test()
{
  cout << "Total time: " << total_ << " seconds\n";
  if(!dir_.empty())
  {
    popd();
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
