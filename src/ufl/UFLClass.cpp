// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLData.h>
#include <dolfin/ufl/UFLIndex.h>
#include <dolfin/ufl/UFLIntegral.h>

#include <dolfin/log/log.h>

#include <iostream>
#include <sstream>

namespace ufl
{

Object::repr_t const Class::default_repr_ = Object::repr_t("Class");
std::string const Class::default_str_ = "class";

//-----------------------------------------------------------------------------
Class::Class() :
    prefixed_name_("", ""),
    args_repr_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name) :
    prefixed_name_(name, ""),
    cpp_proto_(std::make_pair(name, std::string("")),
               std::vector<Object const *>()),
    args_repr_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& pre, std::string const& pos) :
    prefixed_name_(pre, pos),
    cpp_proto_(std::make_pair(pre, pos), std::vector<Object const *>()),
    args_repr_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name, repr_t const& repr) :
    prefixed_name_(name, ""),
    cpp_proto_(make_proto(repr)),
    args_repr_(make_args_repr(repr))
{
  if (prefixed_name_.first != cpp_proto_.first.first
      || prefixed_name_.second != cpp_proto_.first.second)
  {
    dolfin::error(
        "The representation :\n\t" + repr + "\n"
            "is invalid as it does not match the class " + prefixed_name_.first
            + prefixed_name_.second);
  }
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& pre, std::string const& pos, repr_t const& repr) :
    prefixed_name_(pre, pos),
    cpp_proto_(make_proto(repr)),
    args_repr_(make_args_repr(repr))
{
  if (prefixed_name_.first != cpp_proto_.first.first
      || prefixed_name_.second != cpp_proto_.first.second)
  {
    dolfin::error(
        "The representation :\n\t" + repr + "\n"
            "is invalid as it does not match the class " + prefixed_name_.first
            + prefixed_name_.second);
  }
}

//-----------------------------------------------------------------------------
Class::~Class()
= default;

//-----------------------------------------------------------------------------
std::string const Class::name() const
{
  return prefixed_name_.first + prefixed_name_.second;
}

//-----------------------------------------------------------------------------
std::pair<std::string, std::string> const Class::prefixed_name() const
{
  return prefixed_name_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Class::repr() const
{
  return default_repr_;
}

//-----------------------------------------------------------------------------
std::string const& Class::str() const
{
  return default_str_;
}

//-----------------------------------------------------------------------------
std::string Class::make_name(repr_t repr)
{
  std::string str = repr;
  size_t openpos = str.find("(");
  return str.substr(0, openpos);
}

//-----------------------------------------------------------------------------
template<class OBJ>
  Object::repr_t Class::make_repr(std::vector<OBJ const *> const& args) const
  {
    std::stringstream ret;
    if (prefixed_name().second == "") ret << prefixed_name().first << "(";
    else ret << prefixed_name().first;
    
    ret << Object::make_repr(args);
    
    if (prefixed_name().second == "") ret << ")" << prefixed_name().second;
    else ret << prefixed_name().second;
    
    return ret.str();
  }

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  return Class::make_repr(args);
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1,
                                      Object const *& arg2) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  args.push_back(arg2);
  return Class::make_repr(args);
}

//-----------------------------------------------------------------------------
std::vector<Object::repr_t> Class::make_args_repr(repr_t const& repr,
                                                  bool without_pre_pos) const
{
  std::string str = repr;
  std::vector<repr_t> args;
  
  if (str.length() == 0) return args;
  
  size_t star_pos = 0;
  while (star_pos != std::string::npos && star_pos < str.length())
  {
    star_pos = str.find("*");
    if (star_pos != std::string::npos) str.erase(star_pos, 1);
  }
  
  std::vector<char> open_delimiters;
  open_delimiters.push_back('(');
  open_delimiters.push_back('[');
  open_delimiters.push_back('{');
  
  std::vector<char> close_delimiters;
  close_delimiters.push_back(')');
  close_delimiters.push_back(']');
  close_delimiters.push_back('}');
  
  std::vector<size_t> openpositions(3, 0);
  std::vector<size_t> closepositions(3, 0);
  
  std::string pre;
  std::string pos;
  
  //from here, the name (pre and pos) should have been removed
  if (!without_pre_pos) remove_pre_pos(repr, str, pre, pos);
  
  if (str.length() == 0) return args;
  
  for (dolfin::uint i = 0; i < open_delimiters.size(); ++i)
  {
    openpositions[i] = str.find(open_delimiters[i]);
    closepositions[i] = (
        str.rfind(close_delimiters[i]) == std::string::npos ?
            0 : str.rfind(close_delimiters[i]));
  }
  
  const size_t openpos = *std::min_element(openpositions.begin(),
                                           openpositions.end());
  const size_t closepos = *std::max_element(closepositions.begin(),
                                            closepositions.end());
  
  if (str.substr(openpos, 1) == "(" && openpos == 0) args =
      Object::make_args_repr(str.substr(openpos + 1, closepos - openpos - 1),
                             true);
  else args = Object::make_args_repr(str, true);
  
  return args;
}

//-----------------------------------------------------------------------------
Class::CppProto Class::make_proto(Object::repr_t repr) const
{
  std::string str;
  std::string pre;
  std::string pos;
  
  remove_pre_pos(repr, str, pre, pos);
  std::vector<repr_t> args = Class::make_args_repr(str, true);
  return CppProto(std::make_pair(pre, pos), Object::make_args(args));
}

//-----------------------------------------------------------------------------
Object::repr_t const& Class::arg(size_t i)
{
  if (i >= args_repr_.size())
  {
    for (size_t k = 0; k < args_repr_.size(); ++k)
    {
      std::cout << args_repr_[k] << ", ";
    }
    dolfin::error("Invalid index required for argument.");
  }
  return args_repr_[i];
}

//-----------------------------------------------------------------------------
std::vector<Object::repr_t> const& Class::args()
{
  return args_repr_;
}

//-----------------------------------------------------------------------------
void Class::remove_pre_pos(repr_t const& repr, std::string& str,
                           std::string& pre, std::string& pos) const
{
  str = repr;
  pre.clear();
  pos.clear();
  
  std::vector<std::string> open_delimiters;
  open_delimiters.push_back("(");
  open_delimiters.push_back("[");
  open_delimiters.push_back("{");
  std::vector<std::string> close_delimiters;
  close_delimiters.push_back(")");
  close_delimiters.push_back("]");
  close_delimiters.push_back("}");
  
  std::vector<size_t> open_delimiter_pos(3, 0);
  std::vector<size_t> close_delimiter_pos(3, 0);
  
  for (dolfin::uint i = 0; i < open_delimiters.size(); ++i)
  {
    open_delimiter_pos[i] = str.find(open_delimiters[i]);
    close_delimiter_pos[i] = (
        str.rfind(close_delimiters[i]) == std::string::npos ?
            0 : str.rfind(close_delimiters[i]));
  }
  
  const size_t openpos = *std::min_element(open_delimiter_pos.begin(),
                                           open_delimiter_pos.end());
  const size_t closepos = *std::max_element(close_delimiter_pos.begin(),
                                            close_delimiter_pos.end());
  
  if (str.substr(openpos, 1) == "(" && openpos > 0)
  {
    pre = str.substr(0, openpos);
    if (pre.find("*") != std::string::npos) str.erase(1, openpos);
    else str.erase(0, openpos);
  }
  else
  {
    pre = str.substr(0, openpos + 1);
    pos = str.substr(closepos);
    std::string substring = str.substr(openpos + 1, closepos - 1);
    str = substring;
  }
}

//-----------------------------------------------------------------------------
void Class::display() const
{
  std::cout << "Class '" << this->name() << "'" << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//----------------------------INSTANTIATIONS-----------------------------------
template Object::repr_t Class::make_repr(
    std::vector<Data const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<Expression const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<FixedIndex const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<Index const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<IndexBase const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<Integral const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<Measure const *> const&) const;
template Object::repr_t Class::make_repr(
    std::vector<Object const *> const&) const;
}
