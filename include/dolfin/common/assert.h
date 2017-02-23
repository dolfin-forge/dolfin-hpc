// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-02-22
//

#ifndef __DOLFIN_ASSERT_H
#define __DOLFIN_ASSERT_H

#include <dolfin/log/log.h>

namespace dolfin
{

#define dolfin_assert_non_null(P) dolfin_assert(P != NULL)

#define dolfin_assert_null(P)     dolfin_assert(P == NULL)

} /* namespace dolfin */

#endif /* __DOLFIN_ASSERT_H */
