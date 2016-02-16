// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:
//

#ifndef __DOLFIN_SYSTEM_H
#define __DOLFIN_SYSTEM_H

#include <dolfin/common/Array.h>

#include <string>

namespace dolfin
{

///
std::string basename(std::string file);

///
void glob(std::string const& pattern, Array<std::string>& matches);


} /* namespace dolfin */

#endif /* __DOLFIN_SYSTEM_H */
