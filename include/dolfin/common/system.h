// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SYSTEM_H
#define __DOLFIN_SYSTEM_H

#include <dolfin/common/Array.h>

#include <string>

namespace dolfin
{

///
auto filename(std::string const& name, std::string const& ext,
                     uint counter = 0, int width = 0) -> std::string;

auto strcounter(uint counter = 0, int width = 0) -> std::string;

///
auto basename(std::string file) -> std::string;

///
auto dirname(std::string file) -> std::string;

///
auto path(std::string p0, std::string const& p1) -> std::string;
auto path(std::string p0, std::string const& p1, std::string const& p2) -> std::string;

///
void glob(std::string const& pattern, Array<std::string>& matches);

///
void mkdir(std::string const& dir);

///
auto stat(std::string const& file) -> bool;

///
auto getcwd() -> std::string;

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
auto dirstack() -> Array<std::string>&;

///
template<class T>
auto human_readable(T value) -> std::string
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
