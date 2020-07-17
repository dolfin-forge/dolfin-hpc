// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_OFF_FILE_H
#define __DOLFIN_OFF_FILE_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include "GenericFile.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace dolfin
{

/**
 *  @class  OFFFile
 *
 *  @brief  Provides a reader for surface meshes stored with the Geomview Object
 *          File Format [1].
 *
 *  [1] http://people.sc.fsu.edu/~jburkardt/data/off/off.html
 */

class OFFFile : public GenericFile
{

public:

  ///
  OFFFile(std::string const filename);

  ///
  ~OFFFile() override = default;

  /// Input
  void operator>>(Mesh& mesh) override;

private:

  ///
  void get_next_line(std::ifstream& file, std::string& line);

  /// Simple function to get whitespace separated entries on a line
  template<typename T>
    void split_line(std::string const& line, std::vector<T>& vals)
    {
      vals.clear();
      std::istringstream list(line);
      std::string val;
      while (list >> val)
      {
        T v = 0;
        std::stringstream tmp;
        tmp << val;
        tmp >> v;
        vals.push_back(v);
      }
    }

};

}

#endif
