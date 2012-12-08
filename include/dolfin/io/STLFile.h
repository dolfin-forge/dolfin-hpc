// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifndef __STL_FILE_H
#define __STL_FILE_H

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
    void operator>> (Mesh& mesh);    

  };
  
}

#endif
