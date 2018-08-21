// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// Modified by Aurélien Larcher, 2013-2017.
//
// First added:  2005-11-26

#ifndef __DOLFIN_SPACE_TIME_FUNCTION_H
#define __DOLFIN_SPACE_TIME_FUNCTION_H

#include <dolfin/function/Function.h>

namespace dolfin
{

class SpaceTimeFunction : public Function
{

public:

  /// Constructor
  SpaceTimeFunction(std::string const& basename);

  /// Constructor
  SpaceTimeFunction(std::string const& basename, FiniteElementSpace const& space);

  /// Destructor
  ~SpaceTimeFunction();

  /// Load samples
  uint load();

  /// Save sample
  void eval();

  /// Save sample
  void save();

  /// Save sample from external function
  void save(Function const& function);

  ///
  void disp() const;

  //---------------------------------------------------------------------------

  /// Time dependency
  inline Function& operator()(Time const& t)
  {
    return Function::operator ()(t);
  }

  //---------------------------------------------------------------------------

  /// File name format
  static std::string filename(std::string const& basename, uint id);

  //---------------------------------------------------------------------------

private:

  ///
  static void load(real st, std::string const& sname, Function& w);

  ///
  static void save(real st, std::string const& sname, Function& w);

  std::string basename_;
  std::map<real, std::string> samples_;

  //
  Function * W_;
  std::map<real, std::string>::iterator it0_;
  std::map<real, std::string>::iterator it1_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_SPACE_TIME_FUNCTION_H */
