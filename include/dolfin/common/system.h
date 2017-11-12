// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:
//

#ifndef __DOLFIN_SYSTEM_H
#define __DOLFIN_SYSTEM_H

#include <dolfin/common/Array.h>

#include <string>

namespace dolfin
{

///
std::string filename(std::string const& name, std::string const& ext,
                     uint counter = 0, int width = 0);

std::string strcounter(uint counter = 0, int width = 0);

///
std::string basename(std::string file);

///
std::string dirname(std::string file);

///
void glob(std::string const& pattern, Array<std::string>& matches);

///
void mkdir(std::string const& dir);

///
bool stat(std::string const& file);

///
std::string getcwd();

///
void pwd();

///
void cd(std::string const& dir);

///
void mkdircd(std::string const& dir);

///
void pushd(std::string const& dir);

///
void popd();

///
void dirs(int n, std::string const& dir);

///
Array<std::string>& dirstack();

///
template<class T>
std::string human_readable(T value)
{
  // Use IEC units
  static std::string const _unit[] =
    { "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB" };

  uint i = 0;
  std::size_t x(value);
  while (x >> 10)
  {
    ++i;
    x = x >> 10;
  }
  std::stringstream ss;
  ss << x << _unit[i];
  return ss.str();
}

} /* namespace dolfin */

#endif /* __DOLFIN_SYSTEM_H */
