// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
// This file provides values for common constants.

#ifndef __DOLFIN_CONSTANTS_H
#define __DOLFIN_CONSTANTS_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>

using dolfin::real;

real const DOLFIN_EPS        =   3.0e-16;
real const DOLFIN_SQRT_EPS   =   1.0e-8;
real const DOLFIN_PI         =   3.141592653589793238462;

int const DOLFIN_LINELENGTH =   256;
int const DOLFIN_TERM_WIDTH =   80;

#endif
