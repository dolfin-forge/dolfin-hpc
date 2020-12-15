// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_STL_FILE_H
#define __DOLFIN_STL_FILE_H

#include <dolfin/common/types.h>
#include <dolfin/io/GenericFile.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

/**
 * @class STLFile
 * @brief Provides a serial reader for mesh stored in STL binary format.
 *
 *        Reference: http://www.fabbers.com/tech/STL_Format
 *
 */

class STLFile : public GenericFile
{
public:
  STLFile( const std::string filename );

  ~STLFile() override = default;

  // Input
  void operator>>( Mesh & mesh ) override;

  void set_min( Point const & min_ );
  void set_max( Point const & max_ );

private:
  Point min;
  Point max;
};

}

#endif
