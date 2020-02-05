// Copyright (C) 2018 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/log/LogStream.h>

#include <iostream>

namespace dolfin
{

LogStream cout(&std::cout);
LogStream cerr(&std::cerr);
LogStream clog(&std::clog);

} /* namespace dolfin */
