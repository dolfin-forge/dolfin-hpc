// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

#include <dolfin/log/log.h>

#include <iostream>
#include <sstream>

namespace ufl
{

Object::repr_t const Class::default_repr_ = Object::repr_t("Class");
std::string const Class::default_str_ = "class";

//-----------------------------------------------------------------------------
Class::Class() :
    args_repr_(),
    pre_(),
    pos_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name) :
    cpp_proto_(std::make_pair(name, ""), std::vector<Object const *>()),
    args_repr_(),
    pre_(name),
    pos_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& pre, std::string const& pos) :
    cpp_proto_(std::make_pair(pre, pos), std::vector<Object const *>()),
    args_repr_(),
    pre_(pre),
    pos_(pos)
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name, repr_t const& repr) :
    cpp_proto_(make_proto(repr)),
    args_repr_(make_args_repr(repr)),
    pre_(name),
    pos_()
{
  if (pre_ != cpp_proto_.first.first || pos_ != cpp_proto_.first.second)
  {
    dolfin::error("The representation :\n\t" + repr + "\n"
        "is invalid as it does not match the class " + pre_ + pos_);
  }
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& pre, std::string const& pos, repr_t const& repr) :
    cpp_proto_(make_proto(repr)),
    args_repr_(make_args_repr(repr)),
    pre_(pre),
    pos_(pos)
{
//  std::cout << "in Class " << repr << std::endl;
  if (pre_ != cpp_proto_.first.first || pos_ != cpp_proto_.first.second)
  {
    dolfin::error("The representation :\n\t" + repr + "\n"
        "is invalid as it does not match the class " + pre_ + pos_);
  }
}

//-----------------------------------------------------------------------------
Class::~Class()
{
}

//-----------------------------------------------------------------------------
std::string const& Class::name() const
{
  return pre_ + pos_;
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::repr() const
{
  return default_repr_;
}

//-----------------------------------------------------------------------------
std::string const Class::str() const
{
  return default_str_;
}

//-----------------------------------------------------------------------------
std::string const Class::make_name(repr_t repr)
{
  std::string str = repr;
  size_t openpos = str.find("(");
  return str.substr(0, openpos);
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(
    std::vector<Object const *> const& args) const
{
  std::stringstream ret;
  ret << name() << "(";
  ret << Object::make_repr(args);
  ret << ")";
  return ret.str();
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  return Object::make_repr(args);
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1,
                                      Object const *& arg2) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  args.push_back(arg2);
  return Object::make_repr(args);
}

//-----------------------------------------------------------------------------
std::vector<Object::repr_t> const Class::make_args_repr(
    repr_t const& repr, bool const& without_pre_pos) const
{
  std::string str = repr;
//  std::cout << "Class make args repr " << str << std::endl;

  std::vector<char> open_delimiters;
  open_delimiters.push_back('(');
  open_delimiters.push_back('[');
  open_delimiters.push_back('{');

  std::vector<char> close_delimiters;
  close_delimiters.push_back(')');
  close_delimiters.push_back(']');
  close_delimiters.push_back('}');

  std::vector<size_t> openpositions(3,0);
  std::vector<size_t> closepositions(3,0);

  std::string pre;
  std::string pos;

  //from here, the name (pre and pos) should have been removed
  if(!without_pre_pos)
    remove_pre_pos(repr, str, pre, pos);

  for(dolfin::uint i = 0; i<open_delimiters.size(); ++i)
  {
    openpositions[i] = str.find(open_delimiters[i]);  
    closepositions[i] = (str.rfind(close_delimiters[i]) == std::string::npos ? 
        0 : str.rfind(close_delimiters[i]));  
  }

  const size_t openpos = *std::min_element(openpositions.begin(), openpositions.end());
  const size_t closepos = *std::max_element(closepositions.begin(), closepositions.end());
  
  std::vector<repr_t> args;
  if(str.substr(openpos,1) == "(" && openpos == 0) 
    args = Object::make_args_repr(
        str.substr(openpos + 1, closepos - openpos - 1), true);
  else
    args = Object::make_args_repr(
        str, true);

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
void Class::remove_pre_pos(repr_t const& repr, 
    std::string& str, std::string& pre, std::string& pos) const
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

  std::vector<size_t> open_delimiter_pos(3,0);
  std::vector<size_t> close_delimiter_pos(3,0);

  for(dolfin::uint i = 0; i<open_delimiters.size(); ++i)
  {
    open_delimiter_pos[i] = str.find(open_delimiters[i]);
    close_delimiter_pos[i] = (str.rfind(close_delimiters[i]) == std::string::npos ? 
        0 : str.rfind(close_delimiters[i]));  
  }

  const size_t openpos = *std::min_element(open_delimiter_pos.begin(), open_delimiter_pos.end());
  const size_t closepos = *std::max_element(close_delimiter_pos.begin(), close_delimiter_pos.end());

  if(str.substr(openpos,1) == "(" && openpos > 0)
  {
    pre = str.substr(0, openpos);
    str.erase(0, openpos);
  }
  else
  {
    pre = str.substr(0, openpos+1);
    pos = str.substr(closepos);
    std::string substring = str.substr(openpos+1, closepos-1);
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

//-----------------------------------------------------------------------------
ValueArray::ValueArray() :
    std::vector<dolfin::uint>()
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(dolfin::uint const i) :
    std::vector<dolfin::uint>(1, i)
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(dolfin::uint const k, dolfin::uint const i) :
    std::vector<dolfin::uint>(k, i)
{
}

//-----------------------------------------------------------------------------
ValueArray::~ValueArray()
{
}

//-----------------------------------------------------------------------------
std::string const ValueArray::str() const
{
  std::stringstream ss;
  ss << "(";
  for(ValueArray::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    ss << *it << ",";
  }
  ss << ")";
  return ss.str();
}

}
