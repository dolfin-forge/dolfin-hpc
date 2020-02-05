// Copyright (C) 2019 Julian Hornich.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MAYBE_UNUSED_H
#define __DOLFIN_MAYBE_UNUSED_H

// this define may be used to signal that a variable might be unused in the
// current function. From C++17 on, there is the [[maybe_unused]] attribute,
// that should be used instead.
#define MAYBE_UNUSED(x) (void) x;

#endif // __DOLFIN_MAYBE_UNUSED_H
