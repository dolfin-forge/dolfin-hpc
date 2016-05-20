// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.
//
// Imported from licorne

#ifndef DOLFIN_EVOLUTION_TIME_H
#define DOLFIN_EVOLUTION_TIME_H

#include <dolfin/common/types.h>

#include <map>

namespace dolfin
{

class Time
{

public:

  Time(real T0, real T1);

  ~Time();

  /// Return time interval as a pair of real numbers
  std::pair<real, real> const& interval() const;

  /// Return time interval as a pair of real numbers
  int sign() const;

  /// Is the current time in the time interval
  bool is_valid() const;

  /// T0
  real begin() const;

  /// T1
  real end() const;

  /// Measure of time interval
  real measure() const;

  /// Clock
  real const& clock() const;

  /// Elapsed time
  real elapsed() const;

  /// Remaining time
  real remaining() const;

  /// Display current time
  void show() const;

  /// Display basic information
  void disp() const;

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

  std::pair<real, real> const T_;
  int const sign_;

  //
  real t_;

};

} /* namespace dolfin */

#endif /* DOLFIN_EVOLUTION_TIME_H */
