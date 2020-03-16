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

  //--- Insitu ---
  set( "NodeNormal alpha", DOLFIN_PI / 2. );
}

//-----------------------------------------------------------------------------

ParameterSystem::~ParameterSystem()
{
  for ( iterator it = this->begin(); it != this->end(); ++it )
    delete it->second;
}

//-----------------------------------------------------------------------------

Parameter const & ParameterSystem::get(std::string const& key) const
{
  const_iterator p = this->find(key);

  if (p == this->end())
  {
    error("Unknown parameter \"%s\".", key.c_str());
  }

  return *(p->second);
}

//-----------------------------------------------------------------------------

bool ParameterSystem::defined(std::string const& key) const
{
  return (this->count(key) > 0);
}

//-----------------------------------------------------------------------------

void ParameterSystem::disp() const
{
  for (const_iterator it = this->begin(); it != this->end(); ++it)
  {
    std::stringstream ss;
    ss << std::left << std::setw(32) << it->first << " = ";
    switch( it->second->type() )
    {
      case Parameter::bool_t:
        ss << dynamic_cast< parameter<bool> const * >(it->second)->get();
        break;
      case Parameter::int_t:
        ss << dynamic_cast< parameter<int> const * >(it->second)->get();
        break;
      case Parameter::uint_t:
        ss << dynamic_cast< parameter<uint> const * >(it->second)->get();
        break;
      case Parameter::real_t:
        ss << dynamic_cast< parameter<real> const * >(it->second)->get();
        break;
      case Parameter::string_t:
        ss << dynamic_cast< parameter<std::string> const * >(it->second)->get();
        break;
      default:
        ss << "Unknown Parameter";
    }
    message(ss.str());
  }
}

//-----------------------------------------------------------------------------

void ParameterSystem::delete_parameter( std::string const & key )
{
  if ( defined( key ) )
  {
    std::map< std::string, Parameter * >::iterator p = this->find( key );
    delete p->second;
    this->erase( p );
  }
}

//-----------------------------------------------------------------------------
