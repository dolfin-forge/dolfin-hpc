// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2005-12-21

#include <dolfin/parameter/Parameter.h>

#include <dolfin/parameter/ParameterValue.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Parameter::Parameter(int value) :
    value_(new parameter<int>(value)),
    type_(type_int)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(uint value) :
    value_(new parameter<uint>(value)),
    type_(type_uint)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(real value) :
    value_(new parameter<real>(value)),
    type_(type_real)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(bool value) :
    value_(new parameter<bool>(value)),
    type_(type_bool)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(std::string value) :
    value_(new parameter<std::string>(value)),
    type_(type_string)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(const char* value) :
    value_(new parameter<std::string>(std::string(value))),
    type_(type_string)
{
}
//-----------------------------------------------------------------------------
Parameter::Parameter(Parameter const& other) :
    value_(other.value_ ? other.value_->clone() : NULL),
    type_(other.type_)
{
}
//-----------------------------------------------------------------------------
bool Parameter::operator==(Parameter const& other) const
{
  return *(this->value_) == *(other.value_);
}
//-----------------------------------------------------------------------------
bool Parameter::operator!=(Parameter const& other) const
{
  return *(this->value_) != *(other.value_);
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(int value)
{
  *(this->value_) = value;
  return *this;
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(uint value)
{
  *(this->value_) = value;
  return *this;
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(real value)
{
  *(this->value_) = value;
  return *this;
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(bool value)
{
  *(this->value_) = value;
  return *this;
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(std::string value)
{
  *(this->value_) = value;
  return *this;
}
//-----------------------------------------------------------------------------
Parameter const& Parameter::operator=(Parameter const& other)
{
  if (this != &other)
  {
    delete value_;
    value_ = (other.value_ ? other.value_->clone() : NULL);
    type_  = other.type_;
  }
  return *this;
}
//-----------------------------------------------------------------------------
Parameter::~Parameter()
{
  delete value_;
}
//-----------------------------------------------------------------------------
Parameter::operator int() const
{
  return *value_;
}
//-----------------------------------------------------------------------------
Parameter::operator uint() const
{
  return *value_;
}
//-----------------------------------------------------------------------------
Parameter::operator real() const
{
  return *value_;
}
//-----------------------------------------------------------------------------
Parameter::operator bool() const
{
  return *value_;
}
//-----------------------------------------------------------------------------
Parameter::operator std::string() const
{
  return *value_;
}
//-----------------------------------------------------------------------------
Parameter::Type Parameter::type() const
{
  return type_;
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream,
                                      const Parameter& parameter)
{
  switch (parameter.type())
  {
    case Parameter::type_int:
      stream << "[Parameter: value = " << static_cast<int>(parameter)
          << " (int)]";
      break;
    case Parameter::type_uint:
      stream << "[Parameter: value = " << static_cast<uint>(parameter)
          << " (uint)]";
      break;
    case Parameter::type_real:
      stream << "[Parameter: value = " << static_cast<real>(parameter)
          << " (real)]";
      break;
    case Parameter::type_bool:
      if (static_cast<bool>(parameter))
        stream << "[Parameter: value = true (bool)]";
      else
        stream << "[Parameter: value = false (bool)]";
      break;
    case Parameter::type_string:
      stream << "[Parameter: value = \"" << static_cast<std::string>(parameter)
          << "\" (string)]";
      break;
    default:
      error("Unknown parameter type: %d.", parameter.type_);
      break;
  }

  return stream;
}
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                                 const Parameter& parameter)
{
  switch (parameter.type())
  {
    case Parameter::type_int:
      stream << static_cast<int>(parameter);
      break;
    case Parameter::type_uint:
      stream << static_cast<uint>(parameter);
      break;
    case Parameter::type_real:
      stream << static_cast<real>(parameter);
      break;
    case Parameter::type_bool:
      if (static_cast<bool>(parameter))
        stream << true;
      else
        stream << false;
      break;
    case Parameter::type_string:
      stream << static_cast<std::string>(parameter);
      break;
    default:
      error("Unknown parameter type: %d.", parameter.type_);
      break;
  }

  return stream;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
