// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/function/SpaceTimeFunction.h>

#include <dolfin/la/GenericVector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>

#include <sys/stat.h>
#include <unistd.h>

#include <fstream>
#include <iomanip>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction( std::string const & basename )
  : Function()
  , basename_( basename )
  , W_( nullptr )
  , it0_( samples_.end() )
  , it1_( samples_.end() )
{
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction( std::string const &        basename,
                                      FiniteElementSpace const & space )
  : Function( space )
  , basename_( basename )
  , W_( nullptr )
  , it0_( samples_.end() )
  , it1_( samples_.end() )
{
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::~SpaceTimeFunction()
{
  delete[] W_;
}
//-----------------------------------------------------------------------------
auto SpaceTimeFunction::load() -> size_t
{
  if ( this->empty() )
  {
    error( "SpaceTimeFunction : loading sample for uninitialized function" );
  }
  samples_.clear();
  real   st;
  size_t id = 0;
  for ( std::string sname = filename( basename_, id );
        access( sname.c_str(), F_OK ) == 0;
        sname = filename( basename_, ++id ) )
  {
#ifdef ENABLE_MPIIO
    MPI_File            fh;
    MPI::Communicator & comm = this->mesh().topology().comm();
    MPI::file_open( fh, sname, MPI_MODE_RDONLY, comm );
    MPI::file_read_all( fh, st );
    MPI::file_close( fh );
#else
    std::ifstream fp( sname.c_str(), std::ifstream::binary );
    fp.read( ( char * ) &st, sizeof( real ) );
    fp.close();
#endif
    samples_.insert( std::pair< real, std::string >( st, sname ) );
  }
  return samples_.size();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::eval()
{
  if ( W_ == nullptr )
  {
    if ( this->empty() )
    {
      error( "SpaceTimeFunction : cannot evaluate empty function" );
    }
    W_ = new Function[2];
    W_[0].init( this->space() );
    W_[1].init( this->space() );
  }

  if ( samples_.empty() )
  {
    load();
    if ( samples_.empty() )
    {
      error( "SpaceTimeFunction : no sample found" );
    }
  }

  real const                                  t = this->clock();
  _ordered_map< real, std::string >::iterator it1;
  it1 = samples_.upper_bound( t );
  if ( it1 == samples_.end() )
  {
    --it1;
  }
  if ( it1 == samples_.begin() )
  {
    ++it1;
  }
  _ordered_map< real, std::string >::iterator it0;
  it0 = it1;
  --it0;

  // Do not reload twice
  if ( it0 == it1_ || it1 == it0_ )
  {
    W_[0].swap( W_[1] );
    std::swap( it0_, it1_ );
  }
  if ( it0 != it0_ )
  {
    load( it0->first, it0->second, W_[0] );
    it0_ = it0;
  }
  if ( it1 != it1_ )
  {
    load( it1->first, it1->second, W_[1] );
    it1_ = it1;
  }

  // Interpolate
  real const st = this->clock();
  real const t0 = it0->first;
  real const t1 = it1->first;
  /// @todo Round-off errors may violate maximum principles
  real const w0 = ( t1 - st ) / ( t1 - t0 );
  real const w1 = ( st - t0 ) / ( t1 - t0 );
  this->zero();
  this->axpy( w0, W_[0] );
  this->axpy( w1, W_[1] );
  this->sync();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::save()
{
  save( *this );
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::save( Function const & function )
{
  std::string sname( filename( basename_, samples_.size() ) );
  samples_.insert( std::pair< real, std::string >( function.clock(), sname ) );
  if ( function.empty() )
  {
    error( "SpaceTimeFunction : cannot write empty function" );
  }
  save( function.clock(), sname, const_cast< Function & >( function ) );
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::disp() const
{
  section( "SpaceTimeFunction" );
  for ( std::pair< real const, std::string > const & sample : samples_ )
  {
    message( "%.16e : %s", sample.first, sample.second.c_str() );
  }
  end();
}
//-----------------------------------------------------------------------------
auto SpaceTimeFunction::filename( std::string const & basename, size_t id )
  -> std::string
{
  std::ostringstream filename;
  filename << basename << std::setw( 6 ) << std::setfill( '0' ) << id;
#ifdef ENABLE_MPIIO
  filename << ".bin";
#else
  filename << "_" << MPI::rank() << ".bin";
#endif
  return filename.str();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::load( real t, std::string const & sname, Function & w )
{
  real                st;
  size_t              sp;
  size_t              offset[3] = { 0, 0, 0 };
  std::vector< real > values( w.vector().local_size() );

#ifdef ENABLE_MPIIO
  MPI_File            fh;
  MPI_Offset          byte_offset = 0;
  MPI::Communicator & comm        = w.mesh().topology().comm();
  size_t              pe_rank     = w.mesh().topology().comm_rank();
  size_t              pe_size     = w.mesh().topology().comm_size();
  MPI::file_open( fh, sname, MPI_MODE_RDONLY, comm );
  byte_offset += MPI::file_read_all( fh, st );
  byte_offset += MPI::file_read_at_all( fh, sp, byte_offset );
  if ( sp != pe_size )
  {
    error(
      "SpaceTimeFunction : communicator size mismatch %u != %u", sp, pe_size );
  }
  byte_offset +=
    MPI::file_read_at_all( fh,
                           &offset[0],
                           2,
                           byte_offset + pe_rank * 2 * sizeof( size_t ),
                           pe_size * 2 );
  if ( offset[1] != w.vector().local_size() )
  {
    error( "SpaceTimeFunction : local size mismatch" );
  }
  byte_offset += MPI::file_read_at_all( fh, offset[2], byte_offset );
  MPI::file_read_at_all(
    fh, &values[0], offset[1], byte_offset + offset[0] * sizeof( real ) );
  MPI::file_close( fh );
#else
  std::ifstream fp( sname.c_str(), std::ifstream::binary );
  fp.read( ( char * ) &st, sizeof( real ) );
  fp.read( ( char * ) &sp, sizeof( size_t ) );
  fp.read( ( char * ) &offset[0], 3 * sizeof( size_t ) );
  fp.read( ( char * ) &values[0], offset[1] * sizeof( real ) );
  fp.close();
#endif

  if ( t != st )
  {
    error( "SpaceTimeFunction : time stamp mismatch" );
  }

  w.vector().set( &values[0] );
  w.sync();
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::save( real st, std::string const & sname, Function & w )
{
  GenericVector const & x         = w.vector();
  size_t                offset[3] = { x.offset(), x.local_size(), x.size() };
  real *                values    = new real[offset[1]];
  x.get( values );

#ifdef ENABLE_MPIIO
  MPI_File            fh;
  MPI_Offset          byte_offset = 0;
  MPI::Communicator & comm        = w.mesh().topology().comm();
  size_t              pe_rank     = w.mesh().topology().comm_rank();
  size_t              pe_size     = w.mesh().topology().comm_size();
  MPI::file_open( fh, sname, MPI_MODE_WRONLY | MPI_MODE_CREATE, comm );
  if ( pe_rank == 0 )
  {
    MPI::check_error(
      MPI_File_write( fh, &st, sizeof( real ), MPI_BYTE, MPI_STATUS_IGNORE ) );
  }
  byte_offset += sizeof( real );
  if ( pe_rank == 0 )
  {
    MPI::check_error( MPI_File_write_at( fh,
                                         byte_offset,
                                         &pe_size,
                                         sizeof( size_t ),
                                         MPI_BYTE,
                                         MPI_STATUS_IGNORE ) );
  }
  byte_offset += sizeof( size_t );
  byte_offset +=
    MPI::file_write_at_all( fh,
                            &offset[0],
                            2,
                            byte_offset + pe_rank * 2 * sizeof( size_t ),
                            pe_size * 2 );
  if ( pe_rank == 0 )
  {
    MPI::check_error( MPI_File_write_at( fh,
                                         byte_offset,
                                         &offset[2],
                                         sizeof( size_t ),
                                         MPI_BYTE,
                                         MPI_STATUS_IGNORE ) );
  }
  byte_offset += sizeof( size_t );
  MPI::file_write_at_all(
    fh, values, offset[1], byte_offset + offset[0] * sizeof( real ) );
  MPI::file_close( fh );
#else
  size_t sp = 1;
  std::ofstream fp( sname.c_str(), std::ofstream::binary );
  fp.write( ( char * ) &st, sizeof( real ) );
  fp.write( ( char * ) &sp, sizeof( size_t ) );
  fp.write( ( char * ) &offset[0], 3 * sizeof( size_t ) );
  fp.write( ( char * ) values, offset[1] * sizeof( real ) );
  fp.close();
#endif

  delete[] values;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
