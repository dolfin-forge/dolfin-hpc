// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
// This file provides values for common constants.

#pragma once

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>

namespace dolfin
{

constexpr real DOLFIN_EPS      = 3.0e-16;
constexpr real DOLFIN_SQRT_EPS = 1.0e-8;
constexpr real DOLFIN_PI       = 3.141592653589793238462;

constexpr uint DOLFIN_LINELENGTH = 256;
constexpr uint DOLFIN_TERM_WIDTH = 80;

} // namespace dolfin
