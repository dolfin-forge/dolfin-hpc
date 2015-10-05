// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifndef __DOLFIN_STL_FILE_H
#define __DOLFIN_STL_FILE_H

#include <dolfin/common/types.h>
#include "GenericFile.h"

namespace dolfin
{

class STLFile : public GenericFile
{
public:

  STLFile(const std::string filename);
  ~STLFile();

  // Input
  void operator>>(Mesh& mesh);

private:

  struct stl_vertex
  {
    double v[3];
    dolfin::uint index;

    bool operator <(const stl_vertex& other) const
    {
      return ((v[0] < other.v[0])
          || (v[0] == other.v[0]
              && (v[1] < other.v[1] || (v[1] == other.v[1] && v[2] < other.v[2]))));
    }

    bool operator ==(const stl_vertex& other) const
    {
      return (v[0] == other.v[0] && v[1] == other.v[1] && v[2] == other.v[2]);
    }
  };

};

}

#endif
