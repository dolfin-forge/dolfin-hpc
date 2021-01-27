// Copyright (C) 2018 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
#ifndef __DOLFIN_LOG_H
#define __DOLFIN_LOG_H

#include <dolfin/log/LogStream.h>

namespace dolfin
{

//--- Default streams ---------------------------------------------------------

extern LogStream cout;
extern LogStream cerr;
extern LogStream clog;

//-----------------------------------------------------------------------------

void message( std::string msg );

//-----------------------------------------------------------------------------

void message( char const * msg, ... );

//-----------------------------------------------------------------------------

void message( size_t n, std::string msg );

//-----------------------------------------------------------------------------

void message( size_t n, char const * msg, ... );

//-----------------------------------------------------------------------------

void warning( std::string msg );

//-----------------------------------------------------------------------------

void warning( char const * msg, ... );

//-----------------------------------------------------------------------------

void error( std::string msg );

//-----------------------------------------------------------------------------

void error( char const * msg, ... );

//-----------------------------------------------------------------------------

void debug( std::string   file,
            unsigned long line,
            std::string   func,
            std::string   msg );

//-----------------------------------------------------------------------------

void debug( std::string   file,
            unsigned long line,
            std::string   func,
            char const *  msg,
            ... );

//-----------------------------------------------------------------------------

void assertion( std::string   file,
                unsigned long line,
                std::string   func,
                std::string   msg );

//-----------------------------------------------------------------------------

void assertion( std::string   file,
                unsigned long line,
                std::string   func,
                char const *  msg,
                ... );

//-----------------------------------------------------------------------------

void header( char const * msg, ... );

//-----------------------------------------------------------------------------

void header( std::string msg );

//-----------------------------------------------------------------------------

void begin( char const * msg, ... );

//-----------------------------------------------------------------------------

void begin( std::string msg );

//-----------------------------------------------------------------------------

void section( char const * msg, ... );

//-----------------------------------------------------------------------------

void section( std::string msg );

//-----------------------------------------------------------------------------

template < typename T >
inline void prm( char const * k, T v )
{
  cout << std::setw( 24 ) << k << ": " << v << "\n";
}

//-----------------------------------------------------------------------------

void end();

//-----------------------------------------------------------------------------

void skip();

//-----------------------------------------------------------------------------

void timing( char const * task, real t );

//-----------------------------------------------------------------------------

void mark( char const * msg = "" );

//-----------------------------------------------------------------------------

auto verbose() -> int;
auto verbose( int n ) -> int;
auto silence() -> int;

//-----------------------------------------------------------------------------

} /* namespace dolfin */

// Debug macros
#ifdef __GNUG__
#define dolfin_debug( msg )                                 \
  do                                                        \
  {                                                         \
    dolfin::debug( __FILE__, __LINE__, __FUNCTION__, msg ); \
  } while ( false )
#else // __FUNCTION__ is a non-standard GNU extension, disable for all other
      // compilers
#define dolfin_debug( msg )
#endif

#endif /* __DOLFIN_LOG_H */
