#ifndef __DOLFIN_LABEL_H
#define __DOLFIN_LABEL_H

#include <dolfin/common/types.h>

#include <string>

namespace dolfin
{

//---------------------------------------------------------------------------
template < class T >
struct Label : public std::pair< T *, std::string >
{
  Label( T & t, std::string s )
    : std::pair< T *, std::string >( &t, s )
  {
  }
};

//---------------------------------------------------------------------------
template < class T >
struct LabelList : public std::vector< Label< T > >
{
  LabelList()
    : std::vector< Label< T > >()
  {
  }

  LabelList( Label< T > & l )
    : std::vector< Label< T > >( 1, l )
  {
  }

  LabelList( Label< T > l )
    : std::vector< Label< T > >( 1, l )
  {
  }

  LabelList( uint n, Label< T > & l )
    : std::vector< Label< T > >( n, l )
  {
  }

  LabelList( uint n, Label< T > l )
    : std::vector< Label< T > >( n, l )
  {
  }
};

//---------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_LABEL_H */
