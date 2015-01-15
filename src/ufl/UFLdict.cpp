// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLData.h>
#include <dolfin/ufl/UFLdict.h>
#include <dolfin/ufl/UFLIndex.h>

namespace ufl
{
//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict(
    std::vector<std::pair<KEY const *, VALUE const *> > const& map) :
    Class("{", "}"),
    map_(map)
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      typename std::vector<std::pair<KEY const *, VALUE const *> >::const_iterator it = map_.begin();
      ssrepr << "{" << it->first->repr() << ": " << it->second->repr();
      ssstr << "{" << it->first->str() << ": " << it->second->str();
      for(++it; it != map_.end(); ++it)
      {
        ssrepr << ", " << it->first->repr() << ": " << it->second->repr();
        ssstr << ", " << it->first->str() << ": " << it->second->str();
      }
      ssrepr << "}";
      ssstr << "}";
      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict(
    std::pair<KEY const *, VALUE const *> const& map) :
    Class("{", "}"),
    map_(std::vector<std::pair<KEY const*, VALUE const *> >(1,map))
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      typename std::vector<std::pair<KEY const *, VALUE const *> >::const_iterator it = map_.begin();
      ssrepr << "{" << it->first->repr() << ": " << it->second->repr();
      ssstr << "{" << it->first->str() << ": " << it->second->str();
      for(++it; it != map_.end(); ++it)
      {
        ssrepr << ", " << it->first->repr() << ": " << it->second->repr();
        ssstr << ", " << it->first->str() << ": " << it->second->str();
      }
      ssrepr << "}";
      ssstr << "}";
      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict(
    KEY const & key, VALUE const & value) :
    Class("{", "}"),
    map_(std::vector<std::pair<KEY const*, VALUE const *> >(1,std::make_pair(&key, &value)))
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      typename std::vector<std::pair<KEY const *, VALUE const *> >::const_iterator it = map_.begin();
      ssrepr << "{" << it->first->repr() << ": " << it->second->repr();
      ssstr << "{" << it->first->str() << ": " << it->second->str();
      for(++it; it != map_.end(); ++it)
      {
        ssrepr << ", " << it->first->repr() << ": " << it->second->repr();
        ssstr << ", " << it->first->str() << ": " << it->second->str();
      }
      ssrepr << "}";
      ssstr << "}";
      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict(dict const& other_dict) :
    Class("{", "}"),
    map_(other_dict.map())
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      typename std::vector<std::pair<KEY const *, VALUE const *> >::const_iterator it = map_.begin();
      ssrepr << "{" << it->first->repr() << ": " << it->second->repr();
      ssstr << "{" << it->first->str() << ": " << it->second->str();
      for(++it; it != map_.end(); ++it)
      {
        ssrepr << ", " << it->first->repr() << ": " << it->second->repr();
        ssstr << ", " << it->first->str() << ": " << it->second->str();
      }
      ssrepr << "}";
      ssstr << "}";
      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict(repr_t const& repr) :
    Class("{", "}", repr),
    map_(fill_map(args()))
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      if(map_.size()>0)
      {
        typename std::vector<std::pair<KEY const *, VALUE const *> >::const_iterator it = map_.begin();
        ssrepr << "{" << it->first->repr() << ": " << it->second->repr();
        ssstr << "{" << it->first->str() << ": " << it->second->str();
        for(++it; it != map_.end(); ++it)
        {
          ssrepr << ", " << it->first->repr() << ": " << it->second->repr();
          ssstr << ", " << it->first->str() << ": " << it->second->str();
        }
        ssrepr << "}";
        ssstr << "}";
      }
      else
      {
        ssrepr << "{}";
        ssstr << "{}";
      }

      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::dict() :
    Class("{","}"),
    map_()
    {
      std::stringstream ssrepr;
      std::stringstream ssstr;
      ssrepr << "{}";
      ssstr << "{}";
      repr_ = ssrepr.str();
      str_ = ssstr.str();
    }

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dict<KEY,VALUE>::~dict()
{
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> std::vector<Class const*> const 
dict<KEY,VALUE>::operands(std::string const& name) const
{
  std::vector<Class const*> expr0;
  for(dolfin::uint i=0; i<size(); ++i)
  {
    std::vector<Class const*> expr1
      = map_[i].first->operands(name);
    expr0.insert(expr0.end(), expr1.begin(), expr1.end());
  }

  if(name == this->name())
    expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> std::vector<std::vector<Class const*> > const 
dict<KEY,VALUE>::level_operands(std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<std::vector<Class const*> > > new_operands(size());

  dolfin::uint s = 0;
  for(dolfin::uint i=0; i<size(); ++i)
  {
    new_operands[i] = map_[i].first->level_operands(operands);
    s = std::max(s, (dolfin::uint)new_operands[i].size());
  }

  std::vector<std::vector<Class const*> > tmp(s+1);
  std::vector<Class const*> obj0;
  obj0.push_back(this);

  tmp[0] = obj0;
  for(dolfin::uint i=0; i<tmp.size()-1; ++i)
  {
    for(dolfin::uint k=0; k<new_operands.size(); ++k)
    {
      if(i<new_operands[k].size())
        for(dolfin::uint j=0; j<new_operands[k][i].size(); ++j)
        {
          tmp[i+1].push_back(new_operands[k][i][j]);
//          std::cout << new_operands[k][i][j]->name() << std::endl;
        }
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> dolfin::uint const dict<KEY,VALUE>::size() const
{
  dolfin::uint size = map_.size();
  return size;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> std::vector<std::pair<KEY const *, VALUE const *> > const& dict<KEY,VALUE>::map() const
{
  return map_;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> Object::repr_t const dict<KEY,VALUE>::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> std::string const dict<KEY,VALUE>::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> void dict<KEY,VALUE>::display() const
{
  std::cout << "dict of " << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//-----------------------------------------------------------------------------
template <class KEY, class VALUE> std::vector<std::pair<KEY const *, VALUE const *> > const dict<KEY,VALUE>::fill_map(std::vector<repr_t> const& reprs)
{
  std::vector<std::pair<KEY const *, VALUE const *> > map(reprs.size());
  for(dolfin::uint i=0; i<reprs.size(); ++i)
  {
    std::string str = reprs[i];
    std::string delimiter = ": ";
    size_t scpos = 0;
    size_t currpos = 0;
    while (scpos != std::string::npos || scpos < str.size())
    {
      scpos = str.find(delimiter, currpos);
      if(scpos == std::string::npos)//empty, so go to next 
        break;
      else
      {
        Object::repr_t first(str.substr(0,scpos));
        Object::repr_t second(str.substr(scpos+delimiter.length()));
//        KEY const * key = KEY::create(first);
//        VALUE const * val = new VALUE(second);
#ifdef __SUNPRO_CC
        std::pair<KEY const *, VALUE const*> new_item = std::make_pair(KEY::create(first), new const VALUE(second));
#else
        std::pair<KEY const *, VALUE const*> new_item = std::make_pair(KEY::create(first), new VALUE(second));
#endif
        map[i] = new_item;
        currpos = scpos + 1;
      }
    }
  }
  return map;
}

//----------------------------INSTANTIATIONS-----------------------------------
template class dict<Data, type<dolfin::uint> >;
template class dict<IndexBase, type<dolfin::uint> >;
//template class dict<FixedIndex, type<dolfin::uint> >;
//template class dict<type<dolfin::uint>, type<dolfin::uint> >;
} /* namespace ufl */
