
#include <dolfin/log/LogStream.h>

#include <iostream>

namespace dolfin
{

LogStream cout(&std::cout);
LogStream cerr(&std::cerr);
LogStream clog(&std::clog);

} /* namespace dolfin */
