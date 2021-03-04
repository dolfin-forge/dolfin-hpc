// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/parameter/ParameterSystem.h>

#include <dolfin/common/constants.h>

using namespace dolfin;

// Initialize the global parameter database
ParameterSystem ParameterSystem::parameters;

//-----------------------------------------------------------------------------

ParameterSystem::ParameterSystem()
{
//--- Linear algebra ---
#ifdef HAVE_PETSC
  set( "linear algebra backend", "PETSc" );
#elif HAVE_JANPACK
  set( "linear algebra backend", "JANPACK" );
#else
  set( "linear algebra backend", "none" );
#endif

  //--- General parameters ---
  set( "solution file name", "solution.pvd" );

  //--- Parameters for input/output ---
  set( "output_format", "vtk" );

  //--- Parameters for Krylov solvers ---
  set( "Krylov relative tolerance", 1e-10 );
  set( "Krylov absolute tolerance", 1e-20 );
  set( "Krylov divergence limit", 1e4 );
  set( "Krylov maximum iterations", 10000 );
  set( "Krylov GMRES restart", 30 );
  set( "Krylov shift nonzero", 0.0 );
  set( "Krylov report", true );
  set( "Krylov keep PC", false );
  set( "Krylov error on nonconvergence", true );

  //--- Parameters for AMG solvers ---
  set( "AMG relative tolerance", 1e-10 );
  set( "AMG absolute tolerance", 1e-20 );
  set( "AMG maximum iterations", 10000 );
  set( "AMG pre-smoothing steps", 5 );
  set( "AMG post-smoothing steps", 5 );
  set( "AMG theta", 0.15 );
  set( "AMG levels", 20 );
  set( "AMG keep levels", false );

  //--- Parameter for direct (LU) solver ---
  set( "LU report", true );

  //--- Parameter for PDE solver ---
  set( "PDE linear solver", "iterative" );

  //--- Mesh partitioning ---
  set( "report edge cut", false ); // ?
  set( "Mesh read in serial", false );
#if HAVE_PARMETIS
  set( "Mesh partitioner", "parmetis" );
#elif HAVE_ZOLTAN
  set( "Mesh partitioner", "zoltan" );
#endif

  //--- Load balancing ---
  set( "Load balancer report", false );
  set( "Load balancer redistribute", true );

  //--- Parameters for intersection detection ---
  set( "GTS Tolerance", 0.0 ); // Tolerance of GTS BB
  // define size for trianlge tolerance ("is the point within this triangle?")
  set( "Geometrical Tolerance Interval", 0.0 );
  set( "Geometrical Tolerance Triangle", 0.0 );
  set( "Geometrical Tolerance Tetrahedron", 0.0 );
  set( "SubDomain Geometrical Tolerance", 1e-6 );
  set( "SubDomain Intersect Boundary", true );

  //--- Mesh smoothing ---
  set( "Mesh smoothing restricted by rmin", true );

  //--- Insitu ---
  set( "VisIt directory", "" );

  //--- Node Normals ---
  set( "NodeNormal alpha", DOLFIN_PI / 2. );
  set( "NodeNormal restricted", false );
  set( "NodeNormal dump types", false );

  //--- Adaptivity ---
  set( "Adaptivity refinements", 1 );
  // can only be "midpoint" or "vertices"
  set( "Adaptivity projection type", "midpoint" );
}

//-----------------------------------------------------------------------------

ParameterSystem::~ParameterSystem()
{
  for ( ParameterSystem::value_type & param : *this )
    if ( param.second != nullptr )
      delete param.second;
}

//-----------------------------------------------------------------------------

auto ParameterSystem::get( std::string const & key ) const -> Parameter const &
{
  const_iterator p = this->find( key );

  if ( p == this->end() )
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
  }

  return *( p->second );
}

//-----------------------------------------------------------------------------

auto ParameterSystem::get_type( std::string const & key ) -> Parameter::Type
{
  dolfin_assert( this->defined( key ) );
  return ( *this )[key]->type();
}

//-----------------------------------------------------------------------------

auto ParameterSystem::defined( std::string const & key ) const -> bool
{
  return ( this->count( key ) > 0 );
}

//-----------------------------------------------------------------------------

auto ParameterSystem::to_json() const -> std::string
{
  std::stringstream ss;
  ss << "{\n";
  bool first_line = true;
  for ( ParameterSystem::value_type const & param : *this )
  {
    if ( not first_line )
      ss << ",\n";

    ss << "\t\"" << param.first << "\": ";
    switch ( param.second->type() )
    {
      case Parameter::bool_t:
        ss << std::boolalpha << this->get< bool >( param.first );
        break;
      case Parameter::int_t:
        ss << this->get< int >( param.first );
        break;
      case Parameter::size_t_t:
        ss << this->get< size_t >( param.first );
        break;
      case Parameter::real_t:
        ss << this->get< real >( param.first );
        break;
      case Parameter::string_t:
        ss << "\"" << this->get< std::string >( param.first ) << "\"";
        break;
      default:
        ss << "\"Unknown Parameter Type\";";
    }

    first_line = false;
  }
  ss << "\n}\n";
  return ss.str();
}

//-----------------------------------------------------------------------------

auto ParameterSystem::serialize() const -> std::string
{
  std::stringstream ss;
  // Parameters are encoded as: "name";type;"value";
  for ( ParameterSystem::value_type const & param : *this )
  {
    ss << "\"" << param.first << "\";" << param.second->type() << ";";
    switch ( param.second->type() )
    {
      case Parameter::bool_t:
        ss << "\"" << std::boolalpha << this->get< bool >( param.first )
           << "\";";
        break;
      case Parameter::int_t:
        ss << "\"" << this->get< int >( param.first ) << "\";";
        break;
      case Parameter::size_t_t:
        ss << "\"" << this->get< size_t >( param.first ) << "\";";
        break;
      case Parameter::real_t:
        ss << "\"" << this->get< real >( param.first ) << "\";";
        break;
      case Parameter::string_t:
        ss << "\"" << this->get< std::string >( param.first ) << "\";";
        break;
      default:
        ss << "\"Unknown Parameter Type\";";
    }
  }
  return ss.str();
}

//-----------------------------------------------------------------------------

void ParameterSystem::deserialize( std::string const & parameters )
{
  std::string::size_type pos = 0;

  while ( pos < parameters.size() )
  {
    dolfin_assert( parameters[pos] == '\"' );

    std::string name(
      parameters, pos + 1, parameters.find( "\";", pos ) - pos - 1 );
    pos += name.size() + 3;

    dolfin_assert( parameters[pos + 1] == ';' );

    std::string type( parameters, pos, parameters.find( ";", pos ) - pos );
    pos += type.size() + 1;

    dolfin_assert( parameters[pos] == '\"' );

    std::string value(
      parameters, pos + 1, parameters.find( "\";", pos ) - pos - 1 );
    pos += value.size() + 3;

    switch ( std::atoi( type.c_str() ) )
    {
      case Parameter::bool_t:
        if ( value == "true" )
          this->set( name, true );
        else if ( value == "false" )
          this->set( name, false );
        else
          warning( "Invalid bool value: %s = %s", name.c_str(), value.c_str() );
        break;
      case Parameter::int_t:
        this->set( name, std::atoi( value.c_str() ) );
        break;
      case Parameter::size_t_t:
        this->set( name, static_cast< size_t >( std::atoi( value.c_str() ) ) );
        break;
      case Parameter::real_t:
        this->set( name, static_cast< real >( std::atof( value.c_str() ) ) );
        break;
      case Parameter::string_t:
        this->set( name, value );
        break;
      default:
        warning( "\"Unknown Parameter Type\";" );
    }
  }
}

//-----------------------------------------------------------------------------

void ParameterSystem::disp() const
{
  for ( ParameterSystem::value_type const & param : *this )
  {
    std::stringstream ss;
    ss << std::left << std::setw( 32 ) << param.first << " = ";
    switch ( param.second->type() )
    {
      case Parameter::bool_t:
        ss << std::boolalpha << this->get< bool >( param.first );
        break;
      case Parameter::int_t:
        ss << this->get< int >( param.first );
        break;
      case Parameter::size_t_t:
        ss << this->get< size_t >( param.first );
        break;
      case Parameter::real_t:
        ss << this->get< real >( param.first );
        break;
      case Parameter::string_t:
        ss << this->get< std::string >( param.first );
        break;
      default:
        ss << "Unknown Parameter Type";
    }
    message( ss.str() );
  }
}

//-----------------------------------------------------------------------------

void ParameterSystem::delete_parameter( std::string const & key )
{
  if ( defined( key ) )
  {
    _ordered_map< std::string, Parameter * >::iterator p = this->find( key );
    delete p->second;
    this->erase( p );
  }
}

//-----------------------------------------------------------------------------
