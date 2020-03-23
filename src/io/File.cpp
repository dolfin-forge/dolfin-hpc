// Copyright (C) 2002-2006 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <dolfin/io/BinaryFile.h>
#include <dolfin/io/File.h>
#include <dolfin/io/GenericFile.h>
#include <dolfin/io/OFFFile.h>
#include <dolfin/io/STLFile.h>
#include <dolfin/io/VTKFile.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
File::File( const std::string & filename )
  : file_( NULL )
{
  uint const len = filename.size();

  // Choose file type base on suffix.
  if ( ( filename.rfind( ".bin" ) != std::string::npos )
       and filename.rfind( ".bin" ) + std::string( ".bin" ).size() == len )
  {
    file_ = new BinaryFile( filename );
  }
  else if ( ( filename.rfind( ".off" ) != std::string::npos )
            and filename.rfind( ".off" ) + std::string( ".off" ).size() == len )
  {
    file_ = new OFFFile( filename );
  }
  else if ( ( filename.rfind( ".pvd" ) != std::string::npos )
            and filename.rfind( ".pvd" ) + std::string( ".pvd" ).size() == len )
  {
    file_ = new VTKFile( filename );
  }
  else if ( ( filename.rfind( ".stl" ) != std::string::npos )
            and filename.rfind( ".stl" ) + std::string( ".stl" ).size() == len )
  {
    file_ = new STLFile( filename );
  }
	else
  {
    // Could not deduce file type
    error( "Unknown file type for \"%s\".", filename.c_str() );
  }
}
//-----------------------------------------------------------------------------
File::File( const std::string & filename, Type type )
{
  switch ( type )
  {
    case binary:
      file_ = new BinaryFile( filename );
      break;
    case vtk:
      file_ = new VTKFile( filename );
      break;
    default:
      file_ = NULL;
      error( "Unknown file type for \"%s\".", filename.c_str() );
      break;
  }
}
//-----------------------------------------------------------------------------
File::File( const std::string & filename, real const & t )
  : file_( NULL )
{

  if ( filename.rfind( ".pvd" ) != filename.npos )
  {
    file_ = new VTKFile( filename, t );
  }
  else if ( filename.rfind( ".bin" ) != filename.npos )
  {
    file_ = new BinaryFile( filename, t );
  }
  else
  {
    file_ = NULL;
    error( "Unknown file type for time dependent \"%s\".", filename.c_str() );
  }
}
//-----------------------------------------------------------------------------
File::~File()
{
  delete file_;
}
//-----------------------------------------------------------------------------
void File::operator>>( GenericVector & x )
{
  file_->read();

  *file_ >> x;
}
//-----------------------------------------------------------------------------
void File::operator>>( GenericMatrix & A )
{
  file_->read();

  *file_ >> A;
}
//-----------------------------------------------------------------------------
void File::operator>>( Mesh & mesh )
{
  file_->read();

  *file_ >> mesh;
}
//-----------------------------------------------------------------------------
void File::operator>>( MeshFunction< int > & meshfunction )
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>( MeshFunction< uint > & meshfunction )
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>( MeshFunction< real > & meshfunction )
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>( MeshFunction< bool > & meshfunction )
{
  file_->read();

  *file_ >> meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator>>( Function & f )
{
  file_->read();

  *file_ >> f;
}
//-----------------------------------------------------------------------------
void File::operator>>( ParameterList & parameters )
{
  file_->read();

  *file_ >> parameters;
}
//-----------------------------------------------------------------------------
void File::operator>>( LabelList< Function > & list )
{
  file_->read();

  *file_ >> list;
}
//-----------------------------------------------------------------------------
void File::operator<<( GenericVector & x )
{
  file_->write();

  *file_ << x;
}
//-----------------------------------------------------------------------------
void File::operator<<( GenericMatrix & A )
{
  file_->write();

  *file_ << A;
}
//-----------------------------------------------------------------------------
void File::operator<<( Mesh & mesh )
{
  file_->write();

  *file_ << mesh;
}
//-----------------------------------------------------------------------------
void File::operator<<( MeshFunction< int > & meshfunction )
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<( MeshFunction< uint > & meshfunction )
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<( MeshFunction< real > & meshfunction )
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<( MeshFunction< bool > & meshfunction )
{
  file_->write();

  *file_ << meshfunction;
}
//-----------------------------------------------------------------------------
void File::operator<<( Function & u )
{
  file_->write();

  *file_ << u;
}
//-----------------------------------------------------------------------------
void File::operator<<( ParameterList & parameters )
{
  file_->write();

  *file_ << parameters;
}
//-----------------------------------------------------------------------------
void File::operator<<( LabelList< Function > & list )
{
  file_->write();

  *file_ << list;
}
//-----------------------------------------------------------------------------
void File::set_counter( uint new_value )
{
  file_->set_counter( new_value );
}
//-----------------------------------------------------------------------------
std::string File::basename( std::string file )
{
  size_t beg = file.find_last_of( '/' );
  if ( beg != std::string::npos )
  {
    file.erase( 0, beg + 1 );
  }
  size_t pos = file.find( '.' );
  return file.substr( 0, pos );
}
//-----------------------------------------------------------------------------
std::string File::filename( std::string basename, std::string format )
{
  if ( format == "vtk" )
  {
    basename += ".pvd";
  }
  else if ( format == "binary" )
  {
    basename += ".bin";
  }
  else
  {
    error( "filename : unknown file format '%s'", format.c_str() );
  }
  return basename;
}
//-----------------------------------------------------------------------------
std::string File::filename( std::string basename )
{
  return filename( basename, dolfin_get<std::string>( "output_format" ) );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
