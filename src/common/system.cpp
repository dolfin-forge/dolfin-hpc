// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/system.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <sstream>
#include <iomanip>
#include <climits>
#include <cstdlib>
#include <glob.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
std::string filename(std::string const& name, std::string const& ext,
                     uint counter, int width)
{
  std::stringstream ss;
  ss << name << std::setfill('0') << std::setw(width) << counter << ext;
  return ss.str();
}
//-----------------------------------------------------------------------------
std::string strcounter(uint counter, int width)
{
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(width) << counter;
  return ss.str();
}
//-----------------------------------------------------------------------------
std::string basename(std::string file)
{
  size_t beg = file.find_last_of('/');
  if (beg != std::string::npos)
  {
    file.erase(0, beg + 1);
  }
  size_t pos = file.find('.');
  return file.substr(0, pos);
}
//-----------------------------------------------------------------------------
std::string dirname(std::string file)
{
  size_t beg = file.find_last_of('/');
  return file.substr(0, beg);
}
//-----------------------------------------------------------------------------
std::string path(std::string p0, std::string const& p1)
{
  return p0 + "/" + p1;
}
//-----------------------------------------------------------------------------
std::string path(std::string p0, std::string const& p1, std::string const& p2)
{
  return p0 + "/" + p1 + "/" + p2;
}
//-----------------------------------------------------------------------------
void glob(std::string const& pattern, Array<std::string>& matches)
{
  glob_t match;
  ::glob(pattern.c_str(), GLOB_ERR, nullptr, &match);
  for (unsigned int i = 0; i < match.gl_pathc; ++i)
  {
    matches.push_back(std::string(match.gl_pathv[i]));
  }
  ::globfree(&match);
}
//-----------------------------------------------------------------------------
void mkdir(std::string const& dirpath)
{
  struct stat sb;
  if (dolfin::MPI::rank() == 0 && (::stat(dirpath.c_str(), &sb) < 0)
        && ::mkdir(dirpath.c_str(), S_IRWXU) < 0)
  {
    error("Unable to create directory: '%s'", dirpath.c_str());
  }
#if HAVE_MPI
  MPI::check_error( MPI_Barrier(dolfin::MPI::DOLFIN_COMM) );
#endif
  if (dolfin::MPI::rank() > 0 && (::stat(dirpath.c_str(), &sb) < 0)
        && ::mkdir(dirpath.c_str(), S_IRWXU) < 0)
  {
    error("Unable to create directory: '%s'", dirpath.c_str());
  }
  message( 1, "mkdir: %s", dirpath.c_str() );
}
//-----------------------------------------------------------------------------
bool stat(std::string const& dirpath)
{
  struct stat sb;
  return !(::stat(dirpath.c_str(), &sb) < 0);
}
//-----------------------------------------------------------------------------
std::string getcwd()
{
  char curpath[PATH_MAX] = { 0 };
  if (::getcwd(curpath, sizeof(curpath)) == nullptr)
  {
    message("Cannot get current working directory");
  }
  return std::string(curpath);
}
//-----------------------------------------------------------------------------
void pwd()
{
  message("%s", getcwd().c_str());
}
//-----------------------------------------------------------------------------
void cd(std::string const& dirpath)
{
  if (::chdir(dirpath.c_str()) < 0)
  {
    error("Unable to enter directory: '%s'", dirpath.c_str());
  }
  message( 1, "cd: %s", dirpath.c_str() );
#if HAVE_MPI
  MPI::check_error( MPI_Barrier(dolfin::MPI::DOLFIN_COMM) );
#endif
}
//-----------------------------------------------------------------------------
void mkdircd(std::string const& dirpath)
{
  mkdir(dirpath);
  cd(dirpath);
}
//-----------------------------------------------------------------------------
void pushd(std::string const& dirpath)
{
  dirstack().push_back( getcwd() );
  mkdircd( dirpath );
}
//-----------------------------------------------------------------------------
void popd()
{
  if(dirstack().size() == 0)
  {
    error("Trying to popd with empty dirstack");
  }

  cd(dirstack().back());
  dirstack().pop_back();
}
//-----------------------------------------------------------------------------
void dirs(int n, std::string& dirname)
{
  // Get directory from top of the stack
  if(n >= 0)
  {
    if(dirstack().size() <= uint(n))
    {
      error("Invalid access to directory stack");
    }
    dirname = *(dirstack().end() - n - 1);
  }
  // Get directory from bottom of the stack
  else
  {
    if(dirstack().size() < uint(n))
    {
      error("Invalid access to directory stack");
    }
    dirname = *(dirstack().begin() - n - 1);
  }
}
//-----------------------------------------------------------------------------
Array<std::string>& dirstack()
{
  static Array<std::string> dirstack_;
  return dirstack_;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
