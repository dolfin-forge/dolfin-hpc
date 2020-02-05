// Copyright (C) 2003-2008 Anders Logg and Jim Tilander.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ASSERT_H
#define __DOLFIN_ASSERT_H

#include <dolfin/config/dolfin_config.h>

#include <string>

#if (DEBUG && !(__GNUG__))
#include <cassert>
#endif

namespace dolfin
{

// Helper function for dolfin_assert macro
void __dolfin_assert(std::string file, unsigned long line, std::string func,
                     char const * msg);

} /* namespace dolfin */

// Assertion, only active if DEBUG is defined
#if (DEBUG && __GNUG__)
#define dolfin_assert(check) do { if ( !(check) ) { dolfin::__dolfin_assert(__FILE__, __LINE__, __FUNCTION__, "(" #check ")"); } } while (false)
#elif DEBUG // __FUNCTION__ is a non-standard GNU extension, use C89 assert
#define dolfin_assert(check) assert(check)
#else
#define dolfin_assert(check)
#endif

#endif /* __DOLFIN_ASSERT_H */
