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

#include <dolfin/common/constants.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/io/File.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#ifdef ENABLE_MPIIO
#include <mpi.h>
#endif

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
    fixed_timestep_(true),
    timestep_(k),
    num_intervals_(std::min(N, uint(std::floor(measure_ / k)))),
    evaluated_(false),
    U0(mesh_),
    U1(mesh_),
    u0_t(0.),
    u1_t(0.),
    u0_t_valid_(false),
    u1_t_valid_(false),
    curr_sample_(0)
{
  this->addFiles(U_files_, Ut.name(), interval, num_intervals_);
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(Function& Ut,
                                     std::pair<real, real> interval, uint N) :
    function_(Ut),
    mesh_(Ut.mesh()),
    timespan_(interval),
    measure_(interval.second - interval.first),
    fixed_timestep_(false),
    timestep_(0.0),
    num_intervals_(N),
    evaluated_(false),
    U0(mesh_),
    U1(mesh_),
    u0_t(0.),
    u1_t(0.),
    u0_t_valid_(false),
    u1_t_valid_(false),
    curr_sample_(0)
{
  this->addFiles(U_files_, Ut.name(), interval);
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(SpaceTimeFunction const& f) :
    function_(f.function_),
    mesh_(f.mesh_),
    timespan_(f.timespan_),
    measure_(f.measure_),
    fixed_timestep_(f.fixed_timestep_),
    timestep_(f.timestep_),
    num_intervals_(f.num_intervals_),
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
void SpaceTimeFunction::clear()
{
  U_files_.clear();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::eval(real t)
{
  if (!evaluated_)
  {
    U0.init(function_.space());
    U1.init(function_.space());
    evaluated_ = true;
  }

  if (U_files_.size() < 2)
  {
    error("The minimum number of sample files is two (one interval).");
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

  // If t == 0-, we need to step forward one
  if(it1 == U_files_.begin())
  {
    ++it1;
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

  message("S0: t = %8f ; sample = %s; w0 = %8f", t0, name0.c_str(), w0);
  message("S1: t = %8f ; sample = %s; w1 = %8f", t1, name1.c_str(), w1);

  // Compute interpolated value
  evaluant().vector() = 0.0;
  evaluant().vector().axpy(w0, U0.vector());
  evaluant().vector().axpy(w1, U1.vector());
  evaluant().sync_ghosts();

  dolfin_assert(evaluant().vector().max() == evaluant().vector().max() );
  dolfin_assert(evaluant().vector().min() == evaluant().vector().min() );
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::disp() const
{
  cout << "SpaceTimeFunction" << endl;
  cout << "-----------------" << endl;

  // Begin indentation
  begin("");
  cout << "Name                   : " << function_.name() << endl;
  cout << "Number of intervals    : " << num_intervals_ << endl;
  cout << "Number of sample files : " << (uint) U_files_.size() << endl;
  if(!U_files_.empty())
  {
    begin("List of samples        : ");
    for (std::map<real, std::string>::const_iterator it = U_files_.begin();
        it != U_files_.end(); ++it)
    {
      cout << it->second << " :  t = " << real(it->first) << endl;
    }
    end();
  }
  // End indentation
  end();
  skip();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::write(real t)
{
  if (t + 0.5 * timestep_
      >= measure_ * (real(curr_sample_) / real(num_intervals_)))
  {
    message("Save " + function_.name() + " sample to file");
    std::string filename = getFileName(function_.name(), curr_sample_);
    File binfile(filename);
    binfile << function_.vector();
    addPoint(filename, t);
    ++curr_sample_;
  }
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
void SpaceTimeFunction::getFileList(std::vector<std::string>& filenames,
                                    std::string basename)
{
  std::string curr_filename = getFileName(basename, 0);
  for (uint sample_id = 0; access(curr_filename.c_str(), F_OK) == 0;
      curr_filename = getFileName(basename, ++sample_id))
  {
    filenames.push_back(curr_filename);
  }
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::addFiles(std::map<real, std::string> files,
                                 std::string basename,
                                 std::pair<real, real> interval)
{

#ifdef ENABLE_MPIIO
  std::vector<std::string> filenames;
  this->getFileList(filenames, basename);

  int counter = 0;
  for (std::vector<std::string>::iterator it = filenames.begin();
      it != filenames.end(); ++it)
  {
    std::string filename = *it;

    MPI_File fh;
    MPI_Offset byte_offset;
    BinaryFile::BinaryFileHeader hdr;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFile::BinaryFileHeader), MPI_BYTE,
                      MPI_STATUS_IGNORE);

    byte_offset = sizeof(BinaryFile::BinaryFileHeader);

    uint nfunc;
    MPI_File_read_at_all(fh, byte_offset, &nfunc, sizeof(uint), MPI_BYTE,
                         MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    BinaryFile::BinaryFunctionHeader f_hdr;
    MPI_File_read_at_all(fh, byte_offset, &f_hdr,
                         sizeof(BinaryFile::BinaryFunctionHeader),
                         MPI_BYTE, MPI_STATUS_IGNORE);

    // Temporary load function, and parse time stamp
    addPoint(filename, f_hdr.t);
    ++counter;

    MPI_File_close(&fh);
  }
#else
  error("MPI I/O required for space time functions with arbitrary time step");
#endif

}

//-----------------------------------------------------------------------------
void SpaceTimeFunction::addFiles(std::map<real, std::string> files,
                                 std::string basename,
                                 std::pair<real, real> interval, uint N)
{
  // Fixed time step only !
  std::vector<std::string> filenames;
  this->getFileList(filenames, basename);

  if(filenames.size() == 0)
  {
    return;
  }

  dolfin_assert(N > 0);
  real const subinterval_meas = std::fabs(interval.second - interval.first)
      / real(N);

  uint counter = 0;
  for (std::vector<std::string>::iterator it = filenames.begin();
      it != filenames.end(); ++it)
  {
    real t = subinterval_meas * counter;
    addPoint(*it, t);
    ++counter;
  }
  if (counter > N + 1)
  {
    warning("Number of sample files exceeds the number of sampling intervals.");
  }
  if (counter < N + 1)
  {
    warning("Insufficient number of sample files.");
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

