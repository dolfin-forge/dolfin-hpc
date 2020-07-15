// Copyright (C) 2005-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TIME_DEPENDENT_H
#define __DOLFIN_TIME_DEPENDENT_H

#include <dolfin/log/log.h>
#include <dolfin/evolution/Time.h>

namespace dolfin
{

/// Associates an object with time t

class TimeDependent
{

public:

  /// Constructors with internal clock
  TimeDependent();

  /// Constructors with synchronized clock
  TimeDependent(Time const& time);

  /// Copy constructor
  TimeDependent(TimeDependent const& other);

  /// Destructor
  ~TimeDependent();

  ///
  TimeDependent& swap(TimeDependent& other);

  /// With internal clock the current time is synchronized at each call
  /// With a synchronized clock the current time is always fetched from the
  /// associated time instance
  inline TimeDependent const& operator()(Time const& time) const
  {
    if(t_ == nullptr)
    {
      clock_ = time.clock();
    }
    else if(t_ != &time)
    {
      error("TimeDependent : re-associating with another time instance");
    }
    return *this;
  }

  /// Synchronize with another time dependent object, do not expose details.
  inline TimeDependent const& operator()(TimeDependent const& other) const
  {
    if(t_ == nullptr)
    {
      clock_ = other.clock_;
    }
    else if(t_ != other.t_)
    {
      error("TimeDependent : re-associating with another time instance");
    }
    return *this;
  }

  /// Return the time associated with the instance
  inline real clock() const { return (t_ == nullptr ? clock_ : t_->clock()); }

private:

  // Internal clock
  mutable real clock_;

  // Pointer to the current time
  Time const * t_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_TIME_DEPENDENT_H */
