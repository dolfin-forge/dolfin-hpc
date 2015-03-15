// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2007.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2003-03-13
// Last changed: 2014-05-29

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include <dolfin/common/constants.h>
#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/log/Logger.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Logger::Logger() :
    destination_(terminal),
    debug_level_(0),
    indentation_level_(0),
    logstream_(&std::cout)
{
  // Do nothing
  if (dolfin::MPI::processNumber() > 0)
  {
    this->destination_ = silent;
    std::cout.clear(std::iostream::badbit);
  }
}
//-----------------------------------------------------------------------------
Logger::~Logger()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void Logger::message(std::string msg, int debug_level)
{
  write(debug_level, msg);
}
//-----------------------------------------------------------------------------
void Logger::warning(std::string msg)
{
  std::string s = std::string("*** Warning: ") + msg;
  write(0, s);
}
//-----------------------------------------------------------------------------
void Logger::error(std::string msg)
{
  std::string s = std::string("*** Error: ") + msg;
  write(0, s);
  throw std::runtime_error(s);
}
//-----------------------------------------------------------------------------
void Logger::begin(std::string msg, int debug_level)
{
  // Write a message
  message(msg, debug_level);
  ++indentation_level_;
}
//-----------------------------------------------------------------------------
void Logger::skip()
{
  message("\n", debug_level_);
}
//-----------------------------------------------------------------------------
void Logger::end()
{
  if (indentation_level_ > 0)
  {
    --indentation_level_;
  }
}
//-----------------------------------------------------------------------------
void Logger::progress(std::string title, real p)
{
  int N = DOLFIN_TERM_WIDTH - 15;
  int n = static_cast<int>(p * static_cast<real>(N));

  // Print the title
  std::string s = "| " + title;
  for (uint i = 0; i < (N - title.size() - 1); i++)
    s += " ";
  s += "|";
  write(0, s);

  // Print the progress bar
  s = "|";
  for (int i = 0; i < n; i++)
    s += "=";
  if (n > 0 && n < N)
  {
    s += "|";
    n++;
  }
  for (int i = n; i < N; i++)
  {
    s += "-";
  }
  s += "| ";
  std::stringstream line;
  line << std::setiosflags(std::ios::fixed);
  line << std::setprecision(1);
  line << 100.0 * p;
  s += line.str() + "%";
  write(0, s);
}
//-----------------------------------------------------------------------------
void Logger::setOutputDestination(std::string destination)
{

  // Choose output destination
  if (destination == "terminal")
  {
    this->destination_ = terminal;
    std::cout.clear(std::iostream::goodbit);
    logstream_ = &std::cout;
  }
  else if (destination == "silent")
  {
    this->destination_ = silent;
    std::cout.clear(std::iostream::badbit);
    logstream_ = &std::cout;
  }
  else if (destination == "stream")
  {
    warning("Please provide the actual stream. Using terminal instead.");
    this->destination_ = terminal;
    logstream_ = &std::cout;
  }
  else
  {
    this->destination_ = terminal;
    logstream_ = &std::cout;
    message("Unknown output destination, using plain text.");
  }
}
//-----------------------------------------------------------------------------
void Logger::setOutputDestination(std::ostream& ostream)
{
  logstream_ = &ostream;
  this->destination_ = stream;
}
//-----------------------------------------------------------------------------
void Logger::setDebugLevel(int debug_level)
{
  this->debug_level_ = debug_level;
}
//-----------------------------------------------------------------------------
void Logger::timing(std::string task, real elapsed_time)
{
  // Print a message
  std::stringstream line;
  line << "Elapsed time: " << elapsed_time << " (" << task << ")";
  message(line.str(), 1);

  // Store values for summary
  timing_map_t::iterator it = timings_.find(task);
  if (it == timings_.end())
  {
    std::pair<uint, real> timing(1, elapsed_time);
    timings_[task] = timing;
  }
  else
  {
    it->second.first += 1;
    it->second.second += elapsed_time;
  }
}
//-----------------------------------------------------------------------------
void Logger::summary()
{
  if (timings_.size() == 0)
  {
    message("Summary: no timings to report.");
    return;
  }

  message("Summary of timings:");
  for (timing_map_t::iterator it = timings_.begin(); it != timings_.end(); ++it)
  {
    const std::string task = it->first;
    const uint num_timings = it->second.first;
    const real total_time = it->second.second;
    const real average_time = total_time / static_cast<real>(num_timings);

    std::stringstream line;
    line << "  " << task << ": " << total_time << " " << average_time << " "
        << num_timings;
    message(line.str());
  }

  // Clear timings
  timings_.clear();
}
//-----------------------------------------------------------------------------
const std::map<std::string, std::pair<dolfin::uint, dolfin::real> >& Logger::timings() const
{
  return timings_;
}
//-----------------------------------------------------------------------------
void Logger::__debug(std::string msg)
{
  std::string s = std::string("Debug: ") + msg;
  write(0, s);
}
//-----------------------------------------------------------------------------
void Logger::__assert(std::string msg)
{
  std::string s = std::string("*** Assertion ") + msg;
  throw std::runtime_error(s);
}
//-----------------------------------------------------------------------------
void Logger::write(int debug_level, std::string& msg)
{
  // Check debug level
  if (debug_level > this->debug_level_)
    return;

  // Indented;
  *logstream_ << std::setw(2*indentation_level_) << " " << msg << std::endl;
}
//----------------------------------------------------------------------------

}

