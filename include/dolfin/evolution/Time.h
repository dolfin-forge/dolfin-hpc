// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#ifndef DOLFIN_EVOLUTION_TIME_H
#define DOLFIN_EVOLUTION_TIME_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Time
{

public:

  using Interval = std::pair<real, real>;

  Time( real T_start = 0.0, real T_end = 0.0, real T_current = 0.0 );

  Time(Interval I);

  Time(Time const& other) = default;

  ~Time() = default;

  /// Operators
  inline operator real&() { return t_; }
  inline operator real const&() const { return t_; }
  inline Time& operator=(real t) { t_ = t; return *this; }
  inline Time& operator+=(real k) { t_ += k; return *this; }
  inline Time& operator-=(real k) { t_ -= k; return *this; }

  /// Return time interval as a pair of real numbers
  Interval const& interval() const;

  /// Return time sign function
  int sign() const;

  /// Is the current time in the time interval
  bool is_valid(real atol = 0.0) const;

  /// T0
  real begin() const;

  /// T1
  real end() const;

  /// Measure of time interval
  real measure() const;

  /// Clock
  real const& clock() const;

  /// Clock [Fragile]
  real & clock();

  /// Elapsed time
  real elapsed() const;

  /// Remaining time
  real remaining() const;

  /// Elapsed normalized time
  real elapsed_normalized() const;

  /// Remaining normalized time
  real remaining_normalized() const;

  /// Display current time
  void show() const;

  /// Display basic information
  void disp() const;

  /// Step
  inline Time& step(real k ) { t_ +=  sign_ * k; return *this; }

  //--- ITERATOR --------------------------------------------------------------

  class iterator
  {

  public:

    /// Constructor
    iterator()
    {
    }

    /// Destructor
    virtual ~iterator()
    {
    }

    //--- INTERFACE -----------------------------------------------------------

    /// Pre-increment
    virtual Time& operator++() = 0;

    /// Pre-decrement
    virtual Time& operator--() = 0;

    /// Is the iterator valid ?
    virtual bool is_valid() const = 0;

    /// Iteration count
    virtual uint count() const = 0;

    /// Display information
    virtual void disp() const = 0;

    //-------------------------------------------------------------------------

  };

  //---------------------------------------------------------------------------

private:

  Interval const T_;
  int const sign_;

  //
  real t_;

};

} /* namespace dolfin */

#endif /* DOLFIN_EVOLUTION_TIME_H */
