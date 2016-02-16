// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:
//

#include <dolfin/common/system.h>

#include <dolfin/log/log.h>

#include <glob.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
std::string basename(std::string file)
{
  size_t beg = file.find_last_of('/');
  if (beg != std::string::npos)
  {
    file.erase(0, beg + 1);
  }
  size_t pos = file.find('.');
  return file.substr(0, pos);
}
//-----------------------------------------------------------------------------
void glob(std::string const& pattern, Array<std::string>& matches)
{
  glob_t match;
  glob(pattern.c_str(), GLOB_ERR, NULL, &match);
  for (unsigned int i = 0; i < match.gl_pathc; ++i)
  {
    matches.push_back(std::string(match.gl_pathv[i]));
  }
  globfree(&match);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
