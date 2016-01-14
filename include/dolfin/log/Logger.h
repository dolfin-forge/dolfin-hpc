// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Thanks to Jim Tilander for many helpful hints.
//
// Modified by Ola Skavhaug, 2007.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2003-03-13
// Last changed: 2014-05-29

#ifndef __DOLFIN_LOGGER_H
#define __DOLFIN_LOGGER_H

#include <dolfin/common/types.h>

#include <string>
#include <ostream>
#include <map>

#ifdef __APPLE__
#undef __assert
#endif

namespace dolfin
{

class Logger
{
public:

  /// Constructor
  Logger();

  /// Destructor
  ~Logger();

  /// Print message
  void message(std::string msg, int debug_level = 0);

  /// Print warning
  void warning(std::string msg);

  /// Print error message and throw exception
  void error(std::string msg);

  /// Begin task (increase indentation level)
  void begin(std::string msg, int debug_level = 0);

  /// End task (decrease indentation level)
  void end();

  /// Skip line
  void skip();

  /// Draw progress bar
  void progress (std::string title, real p);

  /// Set silent output
  void silence();

  /// Set verbose ouput
  void verbose(uint level = 0);

  /// Set output destination ("terminal" or "silent")
  void setOutputDestination(std::string destination);

  /// Set output destination to stream
  void setOutputDestination(std::ostream& stream);

  /// Set debug level
  void setDebugLevel(int debug_level);

  /// Set debug level
  inline int getDebugLevel() { return debug_level_; }

  //  Set file logging
  void file();

  /// Register timing (for later summary)
  void timing(std::string task, real elapsed_time);

  /// Print summary of timings and tasks, clearing stored timings
  void summary();

  /// Return summary of timings
  const std::map<std::string, std::pair<uint, real> >& timings() const;

  /// Helper function for dolfin_debug macro
  void __debug(std::string msg);

  /// Helper function for dolfin_assert macro
  void __assert(std::string msg);

private:

  typedef std::map<std::string, std::pair<dolfin::uint, real> > timing_map_t;

  // Output destination
  enum Destination {terminal, stream, silent};
  Destination destination_;

  // Write message to current output destination
  void write(int debug_level, std::string& msg);

  // Current debug level
  int debug_level_;

  // Current indentation level
  int indentation_level_;

  // Stream for logging
  std::ostream * logstream_;
  std::ofstream * filestream_;

  // List of timings for tasks, map from string to (num_timings, total_time)
  timing_map_t timings_;

};

}

#endif
