// Copyright (C) 2002-2008 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/GenericFile.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>

#include <fstream>

namespace dolfin
{

//-----------------------------------------------------------------------------

GenericFile::GenericFile( std::string const & type,
                          std::string const & filename )
  : type_( type )
  , filename( filename )
  , opened_read( false )
  , opened_write( false )
  , check_header( false )
  , counter( 0 )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( GenericVector & )
{
  read_not_impl( "Vector" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( GenericMatrix & )
{
  read_not_impl( "Matrix" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( Mesh & )
{
  read_not_impl( "Mesh" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( MeshFunction< int > & )
{
  read_not_impl( "MeshFunction<int>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( MeshFunction< size_t > & )
{
  read_not_impl( "MeshFunction<size_t>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( MeshFunction< real > & )
{
  read_not_impl( "MeshFunction<real>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( MeshFunction< bool > & )
{
  read_not_impl( "MeshFunction<bool>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( Function & )
{
  read_not_impl( "Function" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( ParameterList & )
{
  read_not_impl( "ParameterList" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator>>( LabelList< Function > & )
{
  read_not_impl( "Function" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( GenericVector & )
{
  write_not_impl( "Vector" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( GenericMatrix & )
{
  write_not_impl( "Matrix" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( Mesh & )
{
  write_not_impl( "Mesh" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( MeshFunction< int > & )
{
  write_not_impl( "MeshFunction<int>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( MeshFunction< size_t > & )
{
  write_not_impl( "MeshFunction<size_t>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( MeshFunction< real > & )
{
  write_not_impl( "MeshFunction<real>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( MeshFunction< bool > & )
{
  write_not_impl( "MeshFunction<bool>" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( Function & )
{
  write_not_impl( "Function" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( ParameterList & )
{
  write_not_impl( "ParameterList" );
}

//-----------------------------------------------------------------------------

void GenericFile::operator<<( LabelList< Function > & )
{
  write_not_impl( "Function" );
}

//-----------------------------------------------------------------------------

void GenericFile::read()
{
  opened_read = true;
}

//-----------------------------------------------------------------------------

void GenericFile::write()
{
  if ( !opened_write )
  {
    std::ofstream f( filename.c_str() );
    f.close();
  }
  opened_write = true;
}

//-----------------------------------------------------------------------------

void GenericFile::read_not_impl( const std::string object )
{
  error( "Unable to read objects of type %s from %s files.",
         object.c_str(),
         type_.c_str() );
}

//-----------------------------------------------------------------------------

void GenericFile::write_not_impl( const std::string object )
{
  error( "Unable to write objects of type %s to %s files.",
         object.c_str(),
         type_.c_str() );
}

//-----------------------------------------------------------------------------

void GenericFile::parallel_read_not_impl( const std::string object )
{
  if ( MPI::size() > 1 )
  {
    error( "Unable to read objects of type %s from %s files in parallel.",
           object.c_str(),
           type_.c_str() );
  }
}

//-----------------------------------------------------------------------------

void GenericFile::parallel_write_not_impl( const std::string object )
{
  if ( MPI::size() > 1 )
  {
    error( "Unable to write objects of type %s to %s files in parallel.",
           object.c_str(),
           type_.c_str() );
  }
}

//-----------------------------------------------------------------------------

void GenericFile::set_counter( size_t value )
{
  counter = value;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
