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
  inline auto operator=(real t) -> Time& { t_ = t; return *this; }
  inline auto operator+=(real k) -> Time& { t_ += k; return *this; }
  inline auto operator-=(real k) -> Time& { t_ -= k; return *this; }

  /// Return time interval as a pair of real numbers
  auto interval() const -> Interval const&;

  /// Return time sign function
  auto sign() const -> int;

  /// Is the current time in the time interval
  auto is_valid(real atol = 0.0) const -> bool;

  /// T0
  auto begin() const -> real;

  /// T1
  auto end() const -> real;

  /// Measure of time interval
  auto measure() const -> real;

  /// Clock
  auto clock() const -> real const&;

  /// Clock [Fragile]
  auto clock() -> real &;

  /// Elapsed time
  auto elapsed() const -> real;

  /// Remaining time
  auto remaining() const -> real;

  /// Elapsed normalized time
  auto elapsed_normalized() const -> real;

  /// Remaining normalized time
  auto remaining_normalized() const -> real;

  /// Display current time
  void show( std::string const info = "" ) const;

  /// Display basic information
  void disp() const;

  /// Step
  inline auto step(real k ) -> Time& { t_ +=  sign_ * k; return *this; }

  //--- ITERATOR --------------------------------------------------------------

  class iterator
  {

  public:

    /// Constructor
    iterator() = default;

    /// Destructor
    virtual ~iterator() = default;

    //--- INTERFACE -----------------------------------------------------------

    /// Pre-increment
    virtual auto operator++() -> Time& = 0;

    /// Pre-decrement
    virtual auto operator--() -> Time& = 0;

    /// Is the iterator valid ?
    virtual auto is_valid() const -> bool = 0;

    /// Iteration count
    virtual auto count() const -> uint = 0;

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
