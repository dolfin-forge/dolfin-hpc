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

#include <dolfin/log/log.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/LogManager.h>

#include <stdarg.h>
#include <stdio.h>
#include <signal.h>
#include <iomanip>
#include <sstream>

namespace dolfin
{

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
void message(_msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger().message(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void message(int debug_level, _msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger().message(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
#if __sgi
//-----------------------------------------------------------------------------
void message(std::string msg, ...)
{
  read_str(buffer, msg);
  LogManager::logger().message(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void message(int debug_level, std::string msg, ...)
{
  read_str(buffer, msg);
  LogManager::logger().message(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
#endif
//-----------------------------------------------------------------------------
void warning(std::string msg, ...)
{
#ifndef __sgi
  read(buffer, msg);
#else
  read_str(buffer, msg);
#endif
  LogManager::logger().warning(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void error(std::string msg, ...)
{
#ifndef __sgi
  read(buffer, msg);
#else
  read_str(buffer, msg);
#endif
  LogManager::logger().error(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void begin(_msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger().begin(static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void begin(int debug_level, _msg msg, ...)
{
  read(buffer, msg);
  LogManager::logger().begin(static_cast<std::string>(buffer), debug_level);
}
//-----------------------------------------------------------------------------
void end()
{
  LogManager::logger().end();
}
//-----------------------------------------------------------------------------
void endblock()
{
  LogManager::logger().end();
  LogManager::logger().skip();
}
//-----------------------------------------------------------------------------
void skip()
{
  LogManager::logger().skip();
}
//-----------------------------------------------------------------------------
void header(std::string msg, ...)
{
#ifndef __sgi
  read(buffer, msg);
#else
  read_str(buffer, msg);
#endif
  LogManager::logger().message("**** "+static_cast<std::string>(buffer));
}
//-----------------------------------------------------------------------------
void section(std::string msg, ...)
{
  skip();
  message(msg);
  std::stringstream ss;
  ss << std::setw(msg.size()) << std::setfill('-') << "-" << std::endl;
  begin(ss.str());
}
//-----------------------------------------------------------------------------
void timestamp()
{
  message("[%ul]", time(NULL));
}
//-----------------------------------------------------------------------------
void summary()
{
  LogManager::logger().summary();
}
//-----------------------------------------------------------------------------
const std::map<std::string, std::pair<uint, real> >& timings()
{
  return LogManager::logger().timings();
}
//-----------------------------------------------------------------------------
void __debug(std::string file, unsigned long line,
             std::string function, _msg format, ...)
{
  read(buffer, format);
  std::ostringstream ost;
  ost << file << ":" << line << " in " << function << "()";
  std::string msg = std::string(buffer) + " [at " + ost.str() + "]";
  LogManager::logger().__debug(msg);
}
//-----------------------------------------------------------------------------
void __dolfin_assert(std::string file, unsigned long line,
                     std::string function, _msg format, ...)
{
  read(buffer, format);
  std::ostringstream ost;
  ost << file << ":" << line << " in " << function << "()";
  std::string msg = std::string(buffer) + " [at " + ost.str() + "]";
  LogManager::logger().__assert(msg);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

