// Copyright (C) 2010 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-13
// Last changed: 2014-03-13

#include <dolfin/common/Test.h>

#include <dolfin/common/timing.h>
#include <dolfin/common/system.h>
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
Test::Test(int argc, char *argv[]) :
    Startup(argc, argv),
    btest_(false),
    dir_(""),
    total_(0.0),
    padding_(0)
{
  init(argc, argv);
}

//-----------------------------------------------------------------------------
Test::Test(std::string const& dir) :
    Startup(0, NULL),
    btest_(false),
    dir_(dir),
    total_(0.0),
    padding_(0)
{
  init(0, NULL);
  if(!dir_.empty())
  {
    pushd(dir_);
  }
}

//-----------------------------------------------------------------------------
void Test::init(int argc, char *argv[])
{
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
        if (i < argc && (getopt(argc, argv, ":d:m:b") == -1))
        {
          ++optind;
        }
        break;
      }
  }
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
  padding_ = std::max(padding_, (uint) name.size());
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
  cout << "Elapsed time: " << elapsed_time << " seconds" << endl;
  timings_.back().second = elapsed_time;
  total_ += elapsed_time;
  skip();
}

//-----------------------------------------------------------------------------
Test::~Test()
{
  cout << "Total time: " << total_ << " seconds" << endl;
  if(!dir_.empty())
  {
    popd();
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
