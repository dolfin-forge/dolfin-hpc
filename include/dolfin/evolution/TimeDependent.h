// Copyright (C) 2005-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-09-02
// Last changed: 2008-06-23

#ifndef __DOLFIN_TIME_DEPENDENT_H
#define __DOLFIN_TIME_DEPENDENT_H

#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

class Time;

/// Associates an object with time t

class TimeDependent
{

public:

  /// Constructors
  TimeDependent(Time const& time);

  /// Constructors
  TimeDependent(real const& t);

  /// Destructor
  ~TimeDependent();

  /// Return the current time t
  inline real time() const { return t_; }

private:

  // Pointer to the current time
  real const& t_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_TIME_DEPENDENT_H */
