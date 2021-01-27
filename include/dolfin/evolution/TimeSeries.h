// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#ifndef __DOLFIN_TIME_SERIES_H
#define __DOLFIN_TIME_SERIES_H

#include <dolfin/common/types.h>

#include <fstream>

namespace dolfin
{

class TimeSeries
{

public:
  /// Create time series from data
  TimeSeries( std::string const &     filename,
              std::pair< real, real > interval,
              uint                    N,
              real                    k,
              uint                    degree = 1 );

  /// Destructor
  ~TimeSeries();

  /// Clear list of discrete times, values and reset counter to zero
  void clear();

  /// Evaluate values at time t
  void eval( real t );

  /// Write sample
  void write( real t );

  /// Return array of evaluated values
  auto values() -> real *;

  /// Return number of values
  auto value_size() const -> uint;

  /// Return the number of samples registered
  auto num_samples() const -> uint;

  /// Return the time interval of the data samples
  auto sampling_interval() const -> std::pair< real, real >;

  /// Display basic info
  void disp() const;

private:
  /// Load data from file
  void loadData( std::string const & filename );

  /// Add a space function at time t
  void addPoint( real t );

  std::string const             filename_;
  std::pair< real, real > const timespan_;
  real const                    measure_;
  real const                    timestep_;
  uint const                    degree_;

  // Data attributes
  uint                value_size_;
  std::vector< real > values_;

  // Data interval and sampling
  std::fstream            data_file_;
  std::vector< real >     data_values_;
  std::pair< real, real > data_timespan_;
  uint                    num_intervals_;

  _ordered_map< real, uint > discrete_times_;
  real                       t0_;
  real                       t1_;
  uint                       index_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_TIME_SERIES_H */
