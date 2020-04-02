// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETER_SYSTEM_H
#define __DOLFIN_PARAMETER_SYSTEM_H

#include <dolfin/parameter/Parameter.h>

#include <dolfin/log/log.h>

#include <map>

namespace dolfin
{

/// This class holds a global database of parameters for DOLFIN,
/// implemented as a set of (key, value) pairs. Supported value
/// types are bool, int, unit, real, and string.

class ParameterSystem : private std::map< std::string, Parameter * >
{
public:
  /// Singleton instance of global parameter database
  static ParameterSystem parameters;

  //---------------------------------------------------------------------------
public:
  /// Add parameter
  template < typename T >
  void add( std::string key, T const & value );

  /// Set value of parameter
  template < typename T >
  void set( std::string key, T const & value );

  /// Set value of parameter (char * specialization)
  void set( std::string key, char const * value );

  /// Get value of parameter with given key
  template < typename T >
  T get( std::string const & key ) const;

  Parameter::Type get_type( std::string const & key );

  /// Check if parameter with given key has been defined
  bool defined( std::string const & key ) const;

  // (de)serialize ParameterSystem
  std::string serialize() const;
  void deserialize( std::string const & parameters );

  /// print contents of the ParameterSystem
  void disp() const;

  //---------------------------------------------------------------------------
private:
  /// Constructor
  ParameterSystem();

  ~ParameterSystem();

  /// Get Parameter with given key
  Parameter const & get( std::string const & key ) const;

  /// cleanup Parameter
  void delete_parameter( std::string const & key );
};

//-----------------------------------------------------------------------------
// implementation of tempalte functions
//-----------------------------------------------------------------------------

template < typename T >
inline void ParameterSystem::add( std::string key, T const & value )
{
  set( key, value );
}

//-----------------------------------------------------------------------------

template <>
inline void ParameterSystem::set( std::string key, bool const & value )
{
  delete_parameter( key );
  ( *this )[key] = new parameter< bool >( value, Parameter::bool_t );
}

//-----------------------------------------------------------------------------

template <>
inline void ParameterSystem::set( std::string key, int const & value )
{
  delete_parameter( key );
  ( *this )[key] = new parameter< int >( value, Parameter::int_t );
}

//-----------------------------------------------------------------------------

template <>
inline void ParameterSystem::set( std::string key, uint const & value )
{
  delete_parameter( key );
  ( *this )[key] = new parameter< uint >( value, Parameter::uint_t );
}

//-----------------------------------------------------------------------------

template <>
inline void ParameterSystem::set( std::string key, real const & value )
{
  delete_parameter( key );
  ( *this )[key] = new parameter< real >( value, Parameter::real_t );
}

//-----------------------------------------------------------------------------

template <>
inline void ParameterSystem::set( std::string key, std::string const & value )
{
  delete_parameter( key );
  ( *this )[key] = new parameter< std::string >( value, Parameter::string_t );
}

//-----------------------------------------------------------------------------

inline void ParameterSystem::set( std::string key, char const * value )
{
  set( key, std::string( value ) );
}

//-----------------------------------------------------------------------------

template <>
inline bool ParameterSystem::get( std::string const & key ) const
{
  if ( defined( key ) )
  {
    Parameter const & p = get( key );

    switch ( p.type() )
    {
      case Parameter::bool_t:
        return dynamic_cast< parameter< bool > const * >( &p )->get();
        break;
      case Parameter::int_t:
        return dynamic_cast< parameter< int > const * >( &p )->get();
        break;
      case Parameter::uint_t:
        return dynamic_cast< parameter< uint > const * >( &p )->get();
        break;
      case Parameter::real_t:
      case Parameter::string_t:
      default:
        error( "Unknown parameter \"%s\".", key.c_str() );
        return bool(); // dummy
    }
  }
  else
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
    return bool(); // dummy
  }
}

//-----------------------------------------------------------------------------

template <>
inline int ParameterSystem::get( std::string const & key ) const
{
  if ( defined( key ) )
  {
    Parameter const & p = get( key );

    switch ( p.type() )
    {
      case Parameter::bool_t:
        return dynamic_cast< parameter< bool > const * >( &p )->get();
        break;
      case Parameter::int_t:
        return dynamic_cast< parameter< int > const * >( &p )->get();
        break;
      case Parameter::uint_t:
        return dynamic_cast< parameter< uint > const * >( &p )->get();
        break;
      case Parameter::real_t:
        return dynamic_cast< parameter< real > const * >( &p )->get();
        break;
      case Parameter::string_t:
      default:
        error( "Unknown parameter \"%s\".", key.c_str() );
        return int(); // dummy
    }
  }
  else
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
    return int(); // dummy
  }
}

//-----------------------------------------------------------------------------

template <>
inline uint ParameterSystem::get( std::string const & key ) const
{
  if ( defined( key ) )
  {
    Parameter const & p = get( key );

    switch ( p.type() )
    {
      case Parameter::bool_t:
        return dynamic_cast< parameter< bool > const * >( &p )->get();
        break;
      case Parameter::int_t:
        return dynamic_cast< parameter< int > const * >( &p )->get();
        break;
      case Parameter::uint_t:
        return dynamic_cast< parameter< uint > const * >( &p )->get();
        break;
      case Parameter::real_t:
        return dynamic_cast< parameter< real > const * >( &p )->get();
        break;
      case Parameter::string_t:
      default:
        error( "Unknown parameter \"%s\".", key.c_str() );
        return uint(); // dummy
    }
  }
  else
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
    return uint(); // dummy
  }
}

//-----------------------------------------------------------------------------

template <>
inline real ParameterSystem::get( std::string const & key ) const
{
  if ( defined( key ) )
  {
    Parameter const & p = get( key );

    switch ( p.type() )
    {
      case Parameter::bool_t:
        return dynamic_cast< parameter< bool > const * >( &p )->get();
        break;
      case Parameter::int_t:
        return dynamic_cast< parameter< int > const * >( &p )->get();
        break;
      case Parameter::uint_t:
        return dynamic_cast< parameter< uint > const * >( &p )->get();
        break;
      case Parameter::real_t:
        return dynamic_cast< parameter< real > const * >( &p )->get();
        break;
      case Parameter::string_t:
      default:
        error( "Unknown parameter \"%s\".", key.c_str() );
        return real(); // dummy
    }
  }
  else
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
    return real(); // dummy
  }
}

//-----------------------------------------------------------------------------

template <>
inline std::string ParameterSystem::get( std::string const & key ) const
{
  if ( defined( key ) )
  {
    Parameter const & p = get( key );

    switch ( p.type() )
    {
      case Parameter::string_t:
        return dynamic_cast< parameter< std::string > const * >( &p )->get();
        break;
      case Parameter::bool_t:
      case Parameter::int_t:
      case Parameter::uint_t:
      case Parameter::real_t:
      default:
        error( "Unknown parameter \"%s\".", key.c_str() );
        return std::string(); // dummy
    }
  }
  else
  {
    error( "Unknown parameter \"%s\".", key.c_str() );
    return std::string(); // dummy
  }
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif
