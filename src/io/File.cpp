// Copyright (C) 2002-2006 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005, 2006.
// Modified by Haiko Etzel 2005.
// Modified by Magnus Vikstrom 2007.
// Modified by Nuno Lopes 2008
// Modified by Niclas Jansson, 2008-2012.
//
// First added:  2002-11-12
// Last changed: 2012-05-11

#include <dolfin/config/dolfin_config.h>

#include <string>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/io/File.h>
#include <dolfin/io/GenericFile.h>
#include <dolfin/io/XMLFile.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/io/MatlabFile.h>
#include <dolfin/io/OctaveFile.h>
#include <dolfin/io/VTKFile.h>
#include <dolfin/io/RAWFile.h>
#include <dolfin/io/OFFFile.h>
#include <dolfin/io/STLFile.h>
#include <dolfin/io/XYZFile.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
File::File(const std::string& filename) :
    file_(NULL)
{
  // Choose file type base on suffix.

  // FIXME: Use correct funtion to find the suffix; using rfind() makes
  // FIXME: it essential that the suffixes are checked in the correct order.

  if ( filename.rfind(".xml") != filename.npos )
#ifdef HAVE_XML
    file_ = new XMLFile(filename);
#else
    error("DOLFIN is not built with XML support");
#endif
  else if ( filename.rfind(".xml.gz") != filename.npos )
#ifdef HAVE_XML
    file_ = new XMLFile(filename);
#else
    error("DOLFIN is not built with XML support");
#endif
  else if ( filename.rfind(".bin") != filename.npos)
    file_ = new BinaryFile(filename);
  else if ( filename.rfind(".m") != filename.npos )
    file_ = new OctaveFile(filename);
  else if ( filename.rfind(".off") != filename.npos )
    file_ = new OFFFile(filename);
  else if ( filename.rfind(".pvd") != filename.npos )
    file_ = new VTKFile(filename);
  else if ( filename.rfind(".raw") != filename.npos )
    file_ = new RAWFile(filename);
  else if ( filename.rfind(".stl") != filename.npos )
    file_ = new STLFile(filename);
  else if ( filename.rfind(".xyz") != filename.npos )
    file_ = new XYZFile(filename);
  else
  {
    file_ = NULL;
    error("Unknown file type for \"%s\".", filename.c_str());
  }
}
//-----------------------------------------------------------------------------
File::File(const std::string& filename, Type type)
{
  switch ( type ) {
#ifdef HAVE_XML
  case xml:
    file_ = new XMLFile(filename);
    break;
#endif
  case binary:
    file_ = new BinaryFile(filename);
    break;
  case matlab:
    file_ = new MatlabFile(filename);
    break;
  case octave:
    file_ = new OctaveFile(filename);
    break;
  case vtk:
    file_ = new VTKFile(filename);
    break;
  default:
    file_ = NULL;
    error("Unknown file type for \"%s\".", filename.c_str());
    break;
  }
}
//-----------------------------------------------------------------------------
File::File(const std::string& filename, real const& t) :
    file_(NULL)
{

  if ( filename.rfind(".pvd") != filename.npos )
  {
    file_ = new VTKFile(filename, t);
  }
  else if( filename.rfind(".bin") != filename.npos)
  {
    file_ = new BinaryFile(filename , t);
  }
  else
  {
    file_ = NULL;
    error("Unknown file type for time dependent \"%s\".", filename.c_str());
  }
}
//-----------------------------------------------------------------------------
File::~File()
{
  delete file_;
}
//-----------------------------------------------------------------------------
void File::operator>>(GenericVector& x)
{
  file_->read();

  *file_ >> x;
}
//-----------------------------------------------------------------------------
void File::operator>>(GenericMatrix& A)
{
  file_->read();

  *file_ >> A;
}
//-----------------------------------------------------------------------------
void File::operator>>(Mesh& mesh)
{
  file_->read();

  *file_ >> mesh;
}
//-----------------------------------------------------------------------------
void File::operator>>(MeshFunction<int>& meshfunction)
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>(MeshFunction<unsigned int>& meshfunction)
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>(MeshFunction<double>& meshfunction)
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>(MeshFunction<bool>& meshfunction)
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>(Function& f)
{
  file_->read();

  *file_ >> f;
}
//-----------------------------------------------------------------------------
void File::operator>>(ParameterList& parameters)
{
  file_->read();

  *file_ >> parameters;
}
//-----------------------------------------------------------------------------
void File::operator>>(Graph& graph)
{
  file_->read();

  *file_ >> graph;
}
//-----------------------------------------------------------------------------
void File::operator>>(std::vector<std::pair<Function*, std::string> >& f)
{
  file_->read();

  *file_ >> f;
}
//-----------------------------------------------------------------------------
void File::operator<<(GenericVector& x)
{
  file_->write();

  *file_ << x;
}
//-----------------------------------------------------------------------------
void File::operator<<(GenericMatrix& A)
{
  file_->write();

  *file_ << A;
}
//-----------------------------------------------------------------------------
void File::operator<<(Mesh& mesh)
{
  file_->write();

  *file_ << mesh;
}
//-----------------------------------------------------------------------------
void File::operator<<(MeshFunction<int>& meshfunction)
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<(MeshFunction<unsigned int>& meshfunction)
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<(MeshFunction<double>& meshfunction)
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<(MeshFunction<bool>& meshfunction)
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<(Function& u)
{
  file_->write();

  *file_ << u;
}
//-----------------------------------------------------------------------------
void File::operator<<(ParameterList& parameters)
{
  file_->write();

  *file_ << parameters;
}
//-----------------------------------------------------------------------------
void File::operator<<(Graph& graph)
{
  file_->write();

  *file_ << graph;
}
//-----------------------------------------------------------------------------
void File::operator<<(std::vector<std::pair<Function*, std::string> >& f)
{
  file_->write();

  *file_ << f;
}
//-----------------------------------------------------------------------------
void File::set_counter(uint new_value)
{
  file_->set_counter(new_value);
}
//-----------------------------------------------------------------------------
std::string File::basename(std::string file)
{
  size_t beg = file.find_last_of('/');
  if (beg != std::string::npos)
  {
    file.erase(0, beg + 1);
  }
  size_t pos = file.find('.');
  return file.substr(0, pos);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
