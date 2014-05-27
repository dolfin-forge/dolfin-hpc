// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// Modified by Niclas Jansson, 2013.
// Modified by Aurélien Larcher, 2013-2014.
//
// First added:  2005-11-26
// Last changed: 2013-05-17
//

#include <dolfin/adaptivity/SpaceTimeFunction.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/File.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>

#include <cmath>
#include <iomanip>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(Function& Ut,
                                     std::pair<real, real> interval, uint N,
                                     real k) :
    function_(Ut),
    mesh_(Ut.mesh()),
    timespan_(interval),
    measure_(interval.second - interval.first),
    num_intervals_(N),
    fixed_timestep_(true),
    evaluated_(false),
    U0(mesh_),
    U1(mesh_),
    u0_t(0.),
    u1_t(0.),
    u0_t_valid_(false),
    u1_t_valid_(false),
    curr_sample_(0)
{
  std::vector<std::string> filenames;
  this->getFileList(Ut.name(), N, filenames);
  addFiles(filenames, k);
}

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(Function& Ut,
                                     std::pair<real, real> interval, uint N) :
    function_(Ut),
    mesh_(Ut.mesh()),
    timespan_(interval),
    measure_(interval.second - interval.first),
    num_intervals_(N),
    fixed_timestep_(false),
    evaluated_(false),
    U0(mesh_),
    U1(mesh_),
    u0_t(0.),
    u1_t(0.),
    u0_t_valid_(false),
    u1_t_valid_(false),
    curr_sample_(0)
{
  std::vector<std::string> filenames;
  this->getFileList(Ut.name(), N, filenames);
  addFiles(filenames);
}

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(SpaceTimeFunction const& f) :
    function_(f.function_),
    mesh_(f.mesh_),
    timespan_(f.timespan_),
    measure_(f.measure_),
    num_intervals_(f.num_intervals_),
    fixed_timestep_(f.fixed_timestep_),
    evaluated_(f.evaluated_),
    U0(f.function_),
    U1(f.function_),
    u0_t(f.u0_t),
    u1_t(f.u1_t),
    u0_t_valid_(false),
    u1_t_valid_(false),
    curr_sample_(0)
{
}

//-----------------------------------------------------------------------------
SpaceTimeFunction::~SpaceTimeFunction()
{

}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::eval(real t)
{
  if(!evaluated_)
  {
    U0.init(mesh_,function_.signature());
    U1.init(mesh_,function_.signature());
    evaluated_ = true;
  }

  // NOTE: t is the current time in the primal referential t \in [0,primal_Tend]
  // Find element in U_files so that element < t
  std::map<real, std::string>::iterator it1;
  std::map<real, std::string>::iterator it0;

  // Select it1 such that the time t1 is just after t
  it1 = U_files_.upper_bound(t);

  // If t == T, we need to step back one
  if (it1 == U_files_.end())
  {
    --it1;
  }

  it0 = it1;
  --it0;

  real t0 = (*it0).first;
  real t1 = (*it1).first;

  if (t0 != t0 || t1 != t1)
  {
    error("At least one of the iteration times used for interpolation is NaN.");
  }

  std::string name0 = (*it0).second;
  std::string name1 = (*it1).second;

  if (t0 != u0_t || !u0_t_valid_)
  {
    File file0(name0);
    u0_t_valid_ = true;
    u0_t = t0;
    file0 >> U0.vector();
  }

  if (t1 != u1_t || !u1_t_valid_)
  {
    File file1(name1);
    u1_t_valid_ = true;
    u1_t = t1;
    file1 >> U1.vector();
  }

  // Compute weights (linear Lagrange interpolation)
  real w0 = (t1 - t) / (t1 - t0);
  real w1 = (t - t0) / (t1 - t0);

  cout << "S0: t = " << t0 << "; name0 = " << name0 << "; w0 = " << w0 << endl;
  cout << "S1: t = " << t1 << "; name1 = " << name1 << "; w1 = " << w1 << endl;

  // Compute interpolated value
  evaluant().vector() = 0.0;
  evaluant().vector().axpy(w0, U0.vector());
  evaluant().vector().axpy(w1, U1.vector());
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::write(real t)
{
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::write(Array<Function *> functions,
                              std::pair<real, real> interval, uint n, real t)
{

  for (Array<Function *>::iterator it = functions.begin();
      it != functions.end(); ++it)
  {
    Function * f = *it;
    File binfile(getFileName(f->name(), n));
    binfile << f->vector();
  }
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::addPoint(std::string Uname, real t)
{
  U_files_[t] = Uname;
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::addFiles(std::vector<std::string> filenames)
{

#ifdef ENABLE_MPIIO

  int counter = 0;

  for (std::vector<std::string>::iterator it = filenames.begin();
      it != filenames.end(); ++it)
  {
    std::string filename = *it;

    MPI_File fh;
    MPI_Offset byte_offset;
    BinaryFileHeader hdr;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                      MPI_STATUS_IGNORE);

    byte_offset = sizeof(BinaryFileHeader);

    uint nfunc;
    MPI_File_read_at_all(fh, byte_offset, &nfunc, sizeof(uint), MPI_BYTE,
                         MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    BinaryFunctionHeader f_hdr;
    MPI_File_read_at_all(fh, byte_offset, &f_hdr, sizeof(BinaryFunctionHeader),
                         MPI_BYTE, MPI_STATUS_IGNORE);

    // Temporary load function, and parse time stamp
    addPoint(filename, f_hdr.t);

    MPI_File_close(&fh);

    ++counter;
  }
#else
  error("MPI I/O required for space time functions with arbitrary time step");
#endif

}

//-----------------------------------------------------------------------------
void SpaceTimeFunction::addFiles(std::vector<std::string> filenames, real k)
{
  // Fixed time step only !
  int counter = 0;
  int num_files = filenames.size();

  if (num_files == 1)
  {
    error("Number of files is one, divide by zero foreseen.");
  }
  for (std::vector<std::string>::iterator it = filenames.begin();
      it != filenames.end(); ++it)
  {
    std::string filename = *it;

    // OK guys this is *only* valid if we do the right thing i.e:
    // - num_files is the number of samples
    // - T is the measure of the time interval for solving the dual problem
    //	 i.e [sampling_start_time, primal_end_time]
    real t = k * real(counter) / real(num_files - 1);
    addPoint(filename, t);

    ++counter;
  }
  if (counter == 0)
  {
    error("Counter irremediably stayed stuck at zero.");
  }
  for (std::map<real, std::string>::const_iterator it = U_files_.begin();
      it != U_files_.end(); ++it)
  {
    std::cout << std::setw(4) << it->first << " : " << it->second << std::endl;
  }

}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::getFileList(std::string basename, uint N,
                                    std::vector<std::string>& filenames)
{
  filenames.clear();
  // Let us define N as the number of sampling intervals spanning [T0,T1]
  // Therefore the number of files to be loaded is N+1
  // Precondition: N > 0
  if (N == 0)
  {
    error("Trying to interpolate over zero sampling intervals");
  }

  // This loop is merely constructing
  for (uint sample_id = 0; sample_id <= N; ++sample_id)
  {
    filenames.push_back(getFileName(basename, sample_id));
  }
}
//-----------------------------------------------------------------------------
Mesh& SpaceTimeFunction::mesh()
{
  return mesh_;
}
//-----------------------------------------------------------------------------
Function& SpaceTimeFunction::evaluant()
{
  return function_;
}
//-----------------------------------------------------------------------------
std::string SpaceTimeFunction::getFileName(std::string basename, uint sample)
{
  std::stringstream filename;
  filename << basename << std::setw(6) << std::setfill('0') << sample;
#ifdef ENABLE_MPIIO
  filename << ".bin";
#else
  filename << "_" << MPI::processNumber() << ".bin";
#endif
  return filename.str();
}

}

