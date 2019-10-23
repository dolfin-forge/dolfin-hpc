// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#ifndef __DOLFIN_TIME_SERIES_H
#define __DOLFIN_TIME_SERIES_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

#include <fstream>
#include <map>

namespace dolfin
{

class Mesh;

class TimeSeries
{

public:

  /// Create time series from data
  TimeSeries(std::string const& filename, std::pair<real, real> interval,
             uint N, real k, uint degree = 1);

  /// Destructor
  ~TimeSeries();

  /// Clear list of discrete times, values and reset counter to zero
  void clear();

  /// Evaluate values at time t
  void eval(real t);

  /// Write sample
  void write(real t);

  /// Return array of evaluated values
  real * values();

  /// Return number of values
  uint value_size() const;

  /// Return the number of samples registered
  uint num_samples() const;

  /// Return the time interval of the data samples
  std::pair<real, real> sampling_interval() const;

  /// Display basic info
  void disp() const;

private:

  /// Load data from file
  void loadData(std::string const& filename);

  /// Add a space function at time t
  void addPoint(real t);

  std::string const filename_;
  std::pair<real, real> const timespan_;
  real const measure_;
  real const timestep_;
  uint const degree_;

  // Data attributes
  uint value_size_;
  Array<real> values_;

  // Data interval and sampling
  std::fstream data_file_;
  Array<real> data_values_;
  std::pair<real, real> data_timespan_;
  uint num_intervals_;

  std::map<real, uint> discrete_times_;
  real t0_;
  real t1_;
  uint index_;

};

} /* namespace dolfin */

#endif  /* __DOLFIN_TIME_SERIES_H */
