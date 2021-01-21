// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.

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
  ~SpaceTimeFunction() override;

  /// Load samples
  auto load() -> uint;

  /// Save sample
  void eval();

  /// Save sample
  void save();

  /// Save sample from external function
  void save(Function const& function);

  ///
  void disp() const override;

  //---------------------------------------------------------------------------

  /// Time dependency
  inline auto operator()(Time const& t) -> Function&
  {
    return Function::operator ()(t);
  }

  //---------------------------------------------------------------------------

  /// File name format
  static auto filename(std::string const& basename, uint id) -> std::string;

  //---------------------------------------------------------------------------

private:

  ///
  static void load(real st, std::string const& sname, Function& w);

  ///
  static void save(real st, std::string const& sname, Function& w);

  // overload from Function
  void eval(real*, const real*) const override {}

  std::string basename_;
  _ordered_map<real, std::string> samples_;

  //
  Function * W_;
  _ordered_map<real, std::string>::iterator it0_;
  _ordered_map<real, std::string>::iterator it1_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_SPACE_TIME_FUNCTION_H */
