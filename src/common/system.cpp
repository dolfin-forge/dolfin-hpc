// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/system.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <climits>
#include <cstdlib>
#include <glob.h>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

auto filename( std::string const & name,
               std::string const & ext,
               size_t              counter,
               int                 width ) -> std::string
{
  std::stringstream ss;
  ss << name << std::setfill( '0' ) << std::setw( width ) << counter << ext;
  return ss.str();
}

//-----------------------------------------------------------------------------

auto strcounter( size_t counter, int width ) -> std::string
{
  std::stringstream ss;
  ss << std::setfill( '0' ) << std::setw( width ) << counter;
  return ss.str();
}

//-----------------------------------------------------------------------------

auto basename( std::string file ) -> std::string
{
  size_t beg = file.find_last_of( '/' );
  if ( beg != std::string::npos )
  {
    file.erase( 0, beg + 1 );
  }
  size_t pos = file.find( '.' );
  return file.substr( 0, pos );
}

//-----------------------------------------------------------------------------

auto dirname( std::string file ) -> std::string
{
  size_t beg = file.find_last_of( '/' );
  return file.substr( 0, beg );
}

//-----------------------------------------------------------------------------

auto path( std::string p0, std::string const & p1 ) -> std::string
{
  return p0 + "/" + p1;
}

//-----------------------------------------------------------------------------

auto path( std::string p0, std::string const & p1, std::string const & p2 )
  -> std::string
{
  return p0 + "/" + p1 + "/" + p2;
}

//-----------------------------------------------------------------------------

void glob( std::string const & pattern, std::vector< std::string > & matches )
{
  glob_t match;
  ::glob( pattern.c_str(), GLOB_ERR, nullptr, &match );
  for ( unsigned int i = 0; i < match.gl_pathc; ++i )
  {
    matches.push_back( std::string( match.gl_pathv[i] ) );
  }
  ::globfree( &match );
}

//-----------------------------------------------------------------------------

void mkdir( std::string const & dirpath )
{
  struct stat sb;
  if ( dolfin::MPI::rank() == 0 && ( ::stat( dirpath.c_str(), &sb ) < 0 )
       && ::mkdir( dirpath.c_str(), S_IRWXU ) < 0 )
  {
    error( "Unable to create directory: '%s'", dirpath.c_str() );
  }
#if DOLFIN_HAVE_MPI
  MPI::check_error( MPI_Barrier( dolfin::MPI::DOLFIN_COMM ) );
#endif
  if ( dolfin::MPI::rank() > 0 && ( ::stat( dirpath.c_str(), &sb ) < 0 )
       && ::mkdir( dirpath.c_str(), S_IRWXU ) < 0 )
  {
    error( "Unable to create directory: '%s'", dirpath.c_str() );
  }
  message( 1, "mkdir: %s", dirpath.c_str() );
}

//-----------------------------------------------------------------------------

auto stat( std::string const & dirpath ) -> bool
{
  struct stat sb;
  return !( ::stat( dirpath.c_str(), &sb ) < 0 );
}

//-----------------------------------------------------------------------------

auto getcwd() -> std::string
{
  char curpath[PATH_MAX] = { 0 };
  if ( ::getcwd( curpath, sizeof( curpath ) ) == nullptr )
  {
    message( "Cannot get current working directory" );
  }
  return std::string( curpath );
}

//-----------------------------------------------------------------------------

void pwd()
{
  message( "%s", getcwd().c_str() );
}

//-----------------------------------------------------------------------------

void cd( std::string const & dirpath )
{
  if ( ::chdir( dirpath.c_str() ) < 0 )
  {
    error( "Unable to enter directory: '%s'", dirpath.c_str() );
  }

  message( 1, "cd: %s", dirpath.c_str() );

#if DOLFIN_HAVE_MPI
  MPI::check_error( MPI_Barrier( dolfin::MPI::DOLFIN_COMM ) );
#endif
}

//-----------------------------------------------------------------------------

void mkdircd( std::string const & dirpath )
{
  mkdir( dirpath );
  cd( dirpath );
}

//-----------------------------------------------------------------------------

void pushd( std::string const & dirpath )
{
  dirstack().push_back( getcwd() );
  mkdircd( dirpath );
}

//-----------------------------------------------------------------------------

void popd()
{
  if ( dirstack().size() == 0 )
  {
    error( "Trying to popd with empty dirstack" );
  }

  cd( dirstack().back() );
  dirstack().pop_back();
}

//-----------------------------------------------------------------------------

void dirs( int n, std::string & dirname )
{
  // Get directory from top of the stack
  if ( n >= 0 )
  {
    if ( dirstack().size() <= size_t( n ) )
    {
      error( "Invalid access to directory stack" );
    }
    dirname = *( dirstack().end() - n - 1 );
  }
  // Get directory from bottom of the stack
  else
  {
    if ( dirstack().size() < size_t( n ) )
    {
      error( "Invalid access to directory stack" );
    }
    dirname = *( dirstack().begin() - n - 1 );
  }
}

//-----------------------------------------------------------------------------

auto dirstack() -> std::vector< std::string > &
{
  static std::vector< std::string > dirstack_;
  return dirstack_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
