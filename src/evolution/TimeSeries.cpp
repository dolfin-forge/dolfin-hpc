// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#include <dolfin/evolution/TimeSeries.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/main/MPI.h>

#include <sys/stat.h>
#include <unistd.h>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
TimeSeries::TimeSeries(std::string const& filename,
                       std::pair<real, real> interval, uint N, real k,
                       uint degree) :
    filename_(filename),
    timespan_(interval),
    measure_(interval.second - interval.first),
    timestep_(k),
    degree_(degree),
    value_size_(0.0),
    values_(),
    data_file_(),
    data_values_(),
    data_timespan_(timespan_),
    num_intervals_(std::min(N, uint(std::floor(measure_ / k)))),
    discrete_times_(),
    t0_(0.0),
    t1_(0.0),
    index_(0)
{
  loadData(filename_);
}

//-----------------------------------------------------------------------------
TimeSeries::~TimeSeries()
{
  clear();
  data_file_.close();
}

//-----------------------------------------------------------------------------
void TimeSeries::clear()
{
  discrete_times_.clear();
  data_values_.clear();
  values_.clear();
  index_ = 0;
}

//-----------------------------------------------------------------------------
void TimeSeries::eval(real t)
{
  if (!data_values_.empty())
  {
    // Find element in U_files so that element < t
    _ordered_map<real, uint>::iterator it1;
    _ordered_map<real, uint>::iterator it0;

    // Select it1 such that the time t1 is just after t
    it1 = discrete_times_.upper_bound(t);

    // If t == T, we need to step back one
    if (it1 == discrete_times_.end())
    {
      --it1;
    }

    // If t == 0-, we need to step forward one
    if (it1 == discrete_times_.begin())
    {
      ++it1;
    }

    it0 = it1;
    --it0;

    t0_ = it0->first;
    t1_ = it1->first;

    if (t0_ != t0_ || t1_ != t1_)
    {
      error("One of the iteration times used for interpolation is NaN.");
    }

    if (t0_ < data_timespan_.first || t1_ > data_timespan_.second)
    {
      error("Trying to evaluate data outside the defined time span.");
    }

    uint idx0 = it0->second;
    uint idx1 = it1->second;

    // Compute weights (linear Lagrange interpolation)
    real w0 = (t1_ - t) / (t1_ - t0_);
    real w1 = (t - t0_) / (t1_ - t0_);

    message(1, "TimeSeries S0: t = %8f ; sample = %6d; w0 = %8f", t0_, idx0, w0);
    message(1, "TimeSeries S1: t = %8f ; sample = %6d; w1 = %8f", t1_, idx1, w1);
    // Compute interpolated value for each entry
    for (uint i = 0; i < value_size_; ++i)
    {
      values_[i] = w0 * data_values_[idx0 * value_size_ + i]
          + w1 * data_values_[idx1 * value_size_ + i];
    }
  }
  else
  {
    error("Trying to interpolate time series on empty data.");
  }

}

//-----------------------------------------------------------------------------
void TimeSeries::write(real t)
{
  if (MPI::rank() > 0)
  {
    return;
  }

  if(!data_file_.is_open())
  {
    if(access(filename_.c_str(), F_OK) == 0)
    {
      warning("Appending times series data to existing file '%s'",
              filename_.c_str());
    }
    data_file_.open(filename_.c_str(), std::ios_base::out|std::ofstream::app);
  }

  if (t + 0.5 * timestep_ >= measure_ * (real(index_) / real(num_intervals_)))
  {
    message("Save time series at t = %8f", t);
    std::stringstream ss;
    ss << t;
    for (uint i = 0; i < value_size_; ++i)
    {
      ss << "\t" << values_[i];
    }
    data_file_ << ss.str();
    addPoint(t);
  }
}

//-----------------------------------------------------------------------------
real * TimeSeries::values()
{
  return &values_[0];
}

//-----------------------------------------------------------------------------
uint TimeSeries::value_size() const
{
  dolfin_assert(value_size_ == values_.size());
  return value_size_;
}

//-----------------------------------------------------------------------------
uint TimeSeries::num_samples() const
{
  return discrete_times_.size();
}

//-----------------------------------------------------------------------------
std::pair<real, real> TimeSeries::sampling_interval() const
{
  return data_timespan_;
}

//-----------------------------------------------------------------------------
void TimeSeries::disp() const
{
  section("TimeSeries");
  prm("Time interval"     , timespan_);
  prm("Data time span"    , data_timespan_);
  prm("Value size"        , value_size_);
  prm("Number of samples" , data_timespan_);
  prm("Data time span"    , num_samples());
  end();
}

//-----------------------------------------------------------------------------
void TimeSeries::loadData(std::string const&)
{
  Array<real> times;
  if ((MPI::rank() == 0) && (access(filename_.c_str(), F_OK) == 0))
  {
    data_file_.open(filename_.c_str(), std::ios_base::in);
    // Get number of values per line
    std::string line;
    std::getline(data_file_, line);
    std::string buf;
    std::stringstream ss(line);
    Array<std::string> tokens;
    while (ss >> buf)
    {
      tokens.push_back(buf);
    }
    value_size_ = tokens.size() - 1;
    data_values_.clear();
    values_.resize(value_size_);

    // Load data
    data_file_.clear();
    data_file_.seekg(0);
    real time = 0.0;
    real value = 0.0;
    while (data_file_.good())
    {
      data_file_ >> time;
      times.push_back(time);
      for (uint i = 0; i < value_size_; ++i)
      {
        data_file_ >> value;
        data_values_.push_back(value);
      }
    }
    data_file_.close();
  }

  // Exchange data
  if (MPI::size() > 1)
  {
#ifdef HAVE_MPI
    // Value size and number of samples
    uint data_size[2];
    data_size[0] = times.size();
    data_size[1] = value_size_;
    dolfin_assert(value_size_ == values_.size());
    MPI::check_error( MPI_Bcast(&data_size[0], 2, MPI_UNSIGNED, 0,
                                dolfin::MPI::DOLFIN_COMM) );

    // Discrete times
    int tcount = data_size[0];
    times.resize(tcount);
    MPI::check_error( MPI_Bcast(&times[0], tcount, MPI_DOUBLE, 0,
                                dolfin::MPI::DOLFIN_COMM) );

    // Data values
    int vcount = data_size[0] * data_size[1];
    value_size_ = data_size[1];
    values_.resize(value_size_);
    data_values_.resize(vcount);
    MPI::check_error( MPI_Bcast(&data_values_[0], vcount, MPI_DOUBLE, 0,
                                dolfin::MPI::DOLFIN_COMM) );
#endif
  }

  // Update discrete times
  for ( real const & time : times )
  {
    this->addPoint(time);  // add discrete time
  }

  // Update number of intervals
  dolfin_assert(discrete_times_.size() > 1);dolfin_assert(degree_ > 0);
  num_intervals_ = (discrete_times_.size() - 1) / degree_;
}

//-----------------------------------------------------------------------------
void TimeSeries::addPoint(real t)
{
  discrete_times_[t] = index_++;
  if (t < data_timespan_.first)
  {
    data_timespan_.first = t;
  }
  if (t > data_timespan_.second)
  {
    data_timespan_.second = t;
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
