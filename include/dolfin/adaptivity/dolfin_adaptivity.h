#ifndef __DOLFIN_ADAPTIVITY_H
#define __DOLFIN_ADAPTIVITY_H

// DOLFIN adaptivity interface
#include <dolfin/config/dolfin_config.h>
#include <dolfin/adaptivity/AdaptiveRefinement.h>
#ifdef ENABLE_UFL
#include <dolfin/adaptivity/ufc2/AdaptiveRefinementProjectScalar.h>
#include <dolfin/adaptivity/ufc2/AdaptiveRefinementProjectVector.h>
#else
#include <dolfin/adaptivity/ufc1/AdaptiveRefinementProjectScalar.h>
#include <dolfin/adaptivity/ufc1/AdaptiveRefinementProjectVector.h>
#endif
#include <dolfin/adaptivity/SpaceTimeFunction.h>

#endif
