// Copyright (C) 2018 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
#include <dolfin/log/log.h>

#include <ctime>
#include <sstream>
#include <stdexcept>

namespace dolfin
{

#define simple_output(stream, pre, msg, suf) \
  stream << pre << msg << suf;

#define format_output(stream, pre, msg, suf) \
  va_list aptr; \
  va_start(aptr, msg); \
  stream.format(msg, aptr, pre, suf); \
  va_end(aptr);

#define nformat_output(stream, pre, msg, suf, n) \
  va_list aptr; \
  va_start(aptr, msg); \
  stream.format(msg, aptr, pre, suf, &n); \
  va_end(aptr);

//-----------------------------------------------------------------------------
void message(std::string msg)
{
  simple_output(cout, "", msg, "\n");
}

//-----------------------------------------------------------------------------
void message(char const * msg, ...)
{
  format_output(cout, "", msg, "\n");
}

//-----------------------------------------------------------------------------
void message(uint n, std::string msg)
{
  if ((int) n > cout.verbose()) return;
  simple_output(cout, "", msg, "\n");
}

//-----------------------------------------------------------------------------
void message(uint n, char const * msg, ...)
{
  if ((int) n > cout.verbose()) return;
  format_output(cout, "", msg, "\n");
}

//-----------------------------------------------------------------------------
void warning(std::string msg)
{
  simple_output(cout, "Warning: ", msg, "\n");
}

//-----------------------------------------------------------------------------
void warning(char const * msg, ...)
{
  format_output(cout, "Warning: ", msg, "\n");
}

//-----------------------------------------------------------------------------
void error(std::string msg)
{
  simple_output(cerr, "Error  : ", msg, "\n");
  throw std::runtime_error("runtime error");
}

//-----------------------------------------------------------------------------
void error(char const * msg, ...)
{
  format_output(cerr, "Error  : ", msg, "\n");
  throw std::runtime_error("runtime error");
}

//-----------------------------------------------------------------------------
void debug(std::string file, unsigned long line, std::string func, std::string msg)
{
  std::ostringstream os;
  os << " [at " << file << ":" << line << " in " << func << "()]";
  simple_output(cout, "Debug  : ", msg, os.str().c_str());
}

//-----------------------------------------------------------------------------
void debug(std::string file, unsigned long line, std::string func, char const * msg, ...)
{
  std::ostringstream os;
  os << " [at " << file << ":" << line << " in " << func << "()]";
  format_output(cout, "Debug  : ", msg, os.str().c_str());
}

//-----------------------------------------------------------------------------
void assertion(std::string file, unsigned long line, std::string func, std::string msg)
{
  std::ostringstream os;
  os << " [at " << file << ":" << line << " in " << func << "()]";
  simple_output(cerr, "Assert : ", msg, os.str().c_str());
}

//-----------------------------------------------------------------------------
void assertion(std::string file, unsigned long line, std::string func, char const * msg, ...)
{
  std::ostringstream os;
  os << " [at " << file << ":" << line << " in " << func << "()]";
  format_output(cerr, "Assert : ", msg, os.str().c_str());
}

//-----------------------------------------------------------------------------
void header(char const * msg, ...)
{
  format_output(cout, "******** ", msg, "\n");
}

//-----------------------------------------------------------------------------
void header(std::string msg)
{
  simple_output(cerr, "******** ", msg, "\n");
}

//-----------------------------------------------------------------------------
void begin(char const * msg, ...)
{
  format_output(cout, "", msg, "\n"); ++cout;
}

//-----------------------------------------------------------------------------
void begin(std::string msg)
{
  simple_output(cout, "", msg, "\n"); ++cout;
}

//-----------------------------------------------------------------------------
void end()
{
  --cout;
}

//-----------------------------------------------------------------------------
void section(char const * msg, ...)
{
  format_output(cout, "", msg, "\n");
  cout.nputc(8, '-');
  cout << "\n";
  ++cout;
}

//-----------------------------------------------------------------------------
void section(std::string msg)
{
  simple_output(cout, "", msg, "\n");
  cout.nputc(msg.size() - cout.indentwidth(), '-');
  cout << "\n";
  ++cout;
}

//-----------------------------------------------------------------------------
void skip()
{
  cout << "\n";
}

//-----------------------------------------------------------------------------
void timing(char const * task, real t)
{
  cout << "Elapsed time: " << t << " (" << task << ")\n";
}

//-----------------------------------------------------------------------------
void mark(char const * msg)
{
  message("[%9u] %s", std::time(NULL), msg);
}

//-----------------------------------------------------------------------------
int verbose()
{
  return cout.verbose();
}

//-----------------------------------------------------------------------------
int verbose(int n)
{
  return cout.verbose(n);
}

//-----------------------------------------------------------------------------
int silence()
{
  return cout.silence();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
