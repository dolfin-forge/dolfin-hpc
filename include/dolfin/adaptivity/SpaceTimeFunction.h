// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2013-05-17

#ifndef __SPACE_TIME_FUNCTION_H
#define __SPACE_TIME_FUNCTION_H

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
  SpaceTimeFunction(Function& Ut, std::pair<real, real> interval, uint N, real k);

  /// Create space-time function with adaptive time step
  SpaceTimeFunction(Function& Ut, std::pair<real, real> interval, uint N);

  /// Copy constructor
  SpaceTimeFunction(SpaceTimeFunction const& f);

  /// Destructor
  ~SpaceTimeFunction();

  /// Evaluate function at time t, giving result in Ut
  void eval(real t);

  /// Write sample
  void write(real t);

  /// Write function at time t
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

  /// Add a set of functions with arbitrary time steps
  void addFiles(std::vector<std::string> filenames);

  /// Add a set of functions with fixed time step
  void addFiles(std::vector<std::string> filenames, real k);

  ///
  void getFileList(std::string basename, uint N,
                   std::vector<std::string>& files);

  // Evaluant function
  Function& function_;

  // Mesh associated with function
  Mesh& mesh_;

  // Interval and sampling
  std::pair<real, real>  const timespan_;
  real const measure_;
  uint const num_intervals_;
  bool const fixed_timestep_;

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

#ifdef ENABLE_MPIIO

  // File headers for DOLFIN's binary file format, repeated here
  // until this information is available outside of DOLFIN

  enum Binary_data_t
  {
    BINARY_MESH_DATA,
    BINARY_VECTOR_DATA,
    BINARY_FUNCTION_DATA,
    BINARY_MESH_FUNCTION_DATA
  };

  typedef struct
  {
    uint32_t magic;
    uint32_t bendian;
    uint32_t pe_size;
    Binary_data_t type;
  } BinaryFileHeader;

  typedef struct
  {
    uint32_t dim;
    uint32_t size;
    real t;
    char name[256];
  } BinaryFunctionHeader;

#endif
};

}

#endif
