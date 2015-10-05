// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2013-05-17

#ifndef __DOLFIN_SPACE_TIME_FUNCTION_H
#define __DOLFIN_SPACE_TIME_FUNCTION_H

#include <dolfin/function/Function.h>

#include <ufc.h>

#include <stdint.h>

namespace dolfin
{

class Mesh;

class SpaceTimeFunction
{

public:

  /// Create space-time function with fixed time step
  SpaceTimeFunction(Function& Ut, std::pair<real, real> interval, uint N,
                    real k);

  /// Create space-time function with adaptive time step
  SpaceTimeFunction(Function& Ut, std::pair<real, real> interval, uint N);

  /// Copy constructor
  SpaceTimeFunction(SpaceTimeFunction const& f);

  /// Destructor
  ~SpaceTimeFunction();

  /// Clear list of sample files
  void clear();

  /// Evaluate function at time t, giving result in Ut
  void eval(real t);

  /// Write sample
  void write(real t);

  /// List sample files
  void disp() const;

  /// Write function at time t using proper filename formatting
  static void write(Array<Function *> U, std::pair<real, real> interval, uint N,
                    real t);

  /// Return mesh associated with function
  Mesh& mesh();

  /// Return interpolant function
  Function& evaluant();

private:

  /// Add a space function at time t
  void addPoint(std::string Uname, real t);

  /// Filename
  static std::string getFileName(std::string basename, uint sample);

  /// Add a set of functions with fixed time step
  void addFiles(std::map<real, std::string> files, std::string basename,
                std::pair<real, real> interval, uint N);

  /// Add a set of functions with arbitrary time steps
  void addFiles(std::map<real, std::string> files, std::string basename,
                std::pair<real, real> interval);

  void getFileList(std::vector<std::string>& filenames, std::string basename);

  // Evaluant function
  Function& function_;

  // Mesh associated with function
  Mesh& mesh_;

  // Interval and sampling
  std::pair<real, real>  const timespan_;
  real const measure_;
  bool const fixed_timestep_;
  real const timestep_;
  uint const num_intervals_;

  // Space functions defining the current time interval
  bool evaluated_;
  Function U0;
  Function U1;

  real u0_t;
  real u1_t;

  bool u0_t_valid_;
  bool u1_t_valid_;

  std::map<real, std::string> U_files_;
  uint curr_sample_;

};

}

#endif
