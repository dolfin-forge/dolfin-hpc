// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Thanks to Jim Tilander for many helpful hints.
//
// Modified by Ola Skavhaug, 2007.
// Modified by Niclas Jansson, 2009.
//
// First added:  2003-03-13
// Last changed: 2009-05-04

#include <stdarg.h>
#include <stdio.h>
#include <signal.h>
#include <sstream>
#include <dolfin/common/types.h>
#include <dolfin/common/constants.h>
#include <dolfin/log/LogManager.h>
#include <dolfin/log/log.h>

using namespace dolfin;

// Buffers
static char buffer[DOLFIN_LINELENGTH];

#ifdef __sgi
#define read(buffer, msg) \
  va_list aptr; \
  va_start(aptr, msg); \
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg, aptr); \
  va_end(aptr);

#define read_str(buffer, msg) \
  va_list aptr; \
  va_start(aptr, msg); \
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg.c_str(), aptr);	\
  va_end(aptr);
#else
#define read(buffer, msg) \
  va_list aptr; \
  va_start(aptr, msg); \
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg.c_str(), aptr);	\
  va_end(aptr);
#endif
//-----------------------------------------------------------------------------
void dolfin::message(_msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger.message(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void dolfin::message(int debug_level, _msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger.message(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
#if __sgi
//-----------------------------------------------------------------------------
void dolfin::message(std::string msg, ...)
{
  read_str(buffer, msg);
  LogManager::logger.message(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void dolfin::message(int debug_level, std::string msg, ...)
{
  read_str(buffer, msg);
  LogManager::logger.message(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
#endif
//-----------------------------------------------------------------------------
void dolfin::warning(std::string msg, ...)
{
#ifndef __sgi
  read(buffer, msg);
#else
  read_str(buffer, msg);
#endif
  LogManager::logger.warning(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void dolfin::error(std::string msg, ...)
{
#ifndef __sgi
  read(buffer, msg);
#else
  read_str(buffer, msg);
#endif
  LogManager::logger.error(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void dolfin::begin(_msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger.begin(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void dolfin::begin(int debug_level, _msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger.begin(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
void dolfin::end()
{
  LogManager::logger.end();
}
//-----------------------------------------------------------------------------
void dolfin::summary()
{
  LogManager::logger.summary();
}
//-----------------------------------------------------------------------------
const std::map<std::string, std::pair<dolfin::uint, dolfin::real> >& dolfin::timings()
{
  return LogManager::logger.timings();
}
//-----------------------------------------------------------------------------
void dolfin::__debug(std::string file, unsigned long line,
                     std::string function, _msg format, ...)
{
  read(buffer, format);
  std::ostringstream ost;
  ost << file << ":" << line << " in " << function << "()";
  std::string msg = std::string(buffer) + " [at " + ost.str() + "]";
  LogManager::logger.__debug(msg);
}
//-----------------------------------------------------------------------------
void dolfin::__dolfin_assert(std::string file, unsigned long line,
                      std::string function, _msg format, ...)
{
  read(buffer, format);
  std::ostringstream ost;
  ost << file << ":" << line << " in " << function << "()";
  std::string msg = std::string(buffer) + " [at " + ost.str() + "]";
  LogManager::logger.__assert(msg);
}
//-----------------------------------------------------------------------------
