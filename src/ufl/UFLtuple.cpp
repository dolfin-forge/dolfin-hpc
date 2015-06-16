// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#include <dolfin/ufl/UFLData.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>
#include <dolfin/ufl/UFLIntegral.h>
#include <dolfin/ufl/UFLtuple.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::tuple(T const& obj) :
    Class("(", ")"),
    objects_(std::vector<T const *>(1,&obj)),
    repr_(*this, objects_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::tuple(std::vector<T const *> const& objs) :
    Class("(", ")"),
    objects_(objs),
    repr_(*this, objects_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::tuple(tuple<T> const& other_tuple) :
    Class("(", ")"),
    objects_(other_tuple.objects()),
    repr_(*this, objects_),
    str_("")
  {
//    std::cout << "other objs= " << other_tuple.objects().size() << std::endl;
//    std::cout << "objs= " << objects_.size() << std::endl;
  }

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::tuple(repr_t const& repr) :
    Class("(", ")", repr),
    objects_(fill_objects(args())),
    repr_(*this, objects_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::tuple() :
    Class(),
    objects_(),
    repr_(*this, objects_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  template <class T> tuple<T>::~tuple()
  {
  }

//-----------------------------------------------------------------------------
  template <class T> std::vector<Class const*> const tuple<T>::operands(std::string const& name) const
  {
    std::vector<Class const*> obj0;
    for(dolfin::uint i=0; i<objects_.size(); ++i)
    {
      std::vector<Class const*> obj1 = objects_[i]->operands(name);
      obj0.insert(obj0.end(), obj1.begin(), obj1.end());
    }
//    if(name == this->name())
//      obj0.push_back(this);
    return obj0;
  }

//-----------------------------------------------------------------------------
  template <class T> std::vector<std::vector<Class const*> > const tuple<T>::level_operands(
      std::vector<std::vector<Class const*> > const& operands) const
  {
    std::vector<std::vector<std::vector<Class const*> > > new_operands(objects_.size());

    dolfin::uint size = 0;
    for(dolfin::uint i=0; i<objects_.size(); ++i)
    {
      new_operands[i] = objects_[i]->level_operands(operands);
      size = std::max(size, (dolfin::uint)new_operands[i].size());
    }

    std::vector<std::vector<Class const*> > tmp(size+1);
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
//            std::cout << new_operands[k][i][j]->name() << std::endl;
          }
      }
    }

    return tmp;
  }

//-----------------------------------------------------------------------------
  template <class T> dolfin::uint const tuple<T>::size() const
  {
    dolfin::uint size = objects_.size();
    return size;
  }

//-----------------------------------------------------------------------------
  template <class T> std::vector<T const *> const& tuple<T>::operands() const
  {
    return objects_;
  }

//-----------------------------------------------------------------------------
  template <class T> Object::repr_t const& tuple<T>::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  template <class T> std::string const& tuple<T>::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  template <class T> void tuple<T>::display() const
  {
    std::cout << "tuple of " << std::endl;
    std::cout << ".str : " << this->str() << std::endl;
    std::cout << ".repr: " << this->repr() << std::endl;
  }

//-----------------------------------------------------------------------------
  template <class T> std::vector<T const *> const tuple<T>::fill_objects(std::vector<repr_t> const& reprs)
  {
    std::vector<T const *> objs;
    for(dolfin::uint i=0; i<reprs.size(); ++i)
      objs.push_back(T::create(reprs[i]));

    return objs;
  }

//-----------------------------------------------------------------------------
  template <class T> std::vector<T const *> const& tuple<T>::objects() const
  {
    return objects_;
  }

//----------------------------INSTANTIATIONS-----------------------------------
  template class tuple<Index>;
  template class tuple<IndexBase>;
  template class tuple<FixedIndex>;
  template class tuple<Data>;
  template class tuple<Expression>;
  template class tuple<Measure>;
  template class tuple<Integral>;
} /* namespace ufl */
