// Copyright (C) 2003-2008 Anders Logg and Jim Tilander.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2007.
// Modified by Niclas Jansson, 2009-2015.
//
// First added:  2003-03-13
// Last changed: 2015-07-02

#ifndef __DOLFIN_LOG_H
#define __DOLFIN_LOG_H

#include <dolfin/common/types.h>
#include <dolfin/log/LogManager.h>

#include <stdarg.h>
#include <map>
#include <string>

#if (DEBUG && !(__GNUG__))
#include <cassert>
#endif


namespace dolfin
{

  /// The DOLFIN log system provides the following set of functions for
  /// uniform handling of log messages, warnings and errors. In addition,
  /// macros are provided for debug messages and assertions.
  ///
  /// Only messages with a debug level higher than or equal to the global
  /// debug level are printed (the default being zero). The global debug
  /// level may be controlled by
  ///
  ///    set("debug level", debug_level);
  ///
  /// where debug_level is the desired debug level.
  ///
  /// The output destination can be controlled by
  ///
  ///    set("output destination", destination);
  ///
  /// where destination is one of "terminal" (default) or "silent". Setting
  /// the output destination to "silent" means no messages will be printed.

  static Logger& logm = LogManager::logger();

#ifdef __sgi
#define _msg char*
#else
#define _msg std::string
#endif

  /// Print message
  void message(_msg msg, ...);

  /// Print message
  void message(int debug_level, _msg msg, ...);

#if __sgi
  /// Print message
  void message(std::string msg, ...);

  /// Print message
  void message(int debug_level, std::string msg, ...);
#endif
  /// Print warning
  void warning(std::string msg, ...);

  /// Print error message and throw an exception
  void error(std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(int debug_level,_msg msg, ...);

  /// End task (decrease indentation level)
  void end();

  /// End task (decrease indentation level) and skip a line
  void endblock();

  /// Skip line
  void skip();

  /// Write header line
  void header(std::string msg, ...);

  /// Write section opening line
  void section(std::string msg, ...);

  /// Write time stamp marker
  void timestamp();

  /// Print summary of timings and tasks, clearing stored timings
  void summary();

  /// Return summary of timings
  const std::map<std::string, std::pair<dolfin::uint, dolfin::real> >& timings();

  // Helper function for dolfin_debug macro
  void __debug(std::string file, unsigned long line, std::string function, _msg format, ...);

  // Helper function for dolfin_assert macro
  void __dolfin_assert(std::string file, unsigned long line, std::string function, _msg format, ...);

}

// Debug macros (with varying number of arguments)
#ifdef __GNUG__
#define dolfin_debug(msg)              do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg); } while (false)
#define dolfin_debug1(msg, a0)         do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0); } while (false)
#define dolfin_debug2(msg, a0, a1)     do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1); } while (false)
#define dolfin_debug3(msg, a0, a1, a2) do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1, a2); } while (false)
#else  // __FUNCTION__ is a non-standard GNU extension, disable for all other compilers
#define dolfin_debug(msg)
#define dolfin_debug1(msg, a0)
#define dolfin_debug2(msg, a0, a1)
#define dolfin_debug3(msg, a0, a1, a2)
#endif

// Assertion, only active if DEBUG is defined
#if (DEBUG && __GNUG__)
#define dolfin_assert(check) do { if ( !(check) ) { dolfin::__dolfin_assert(__FILE__, __LINE__, __FUNCTION__, "(" #check ")"); } } while (false)
#elif DEBUG // __FUNCTION__ is a non-standard GNU extension, use C89 assert
#define dolfin_assert(check) assert(check)
#else 
#define dolfin_assert(check)
#endif

#endif /* __DOLFIN_LOG_H */
