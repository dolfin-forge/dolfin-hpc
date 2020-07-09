
#include <dolfin/common/assert.h>

#include <sstream>
#include <stdexcept>

namespace dolfin
{

//-----------------------------------------------------------------------------
void __dolfin_assert( std::string   file,
                      unsigned long line,
                      std::string   func,
                      char const *  msg )
{
  std::ostringstream sout;
  sout << msg << " [at " << file << ":" << line << " in " << func << "()]\n";
  throw std::runtime_error( sout.str() );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
