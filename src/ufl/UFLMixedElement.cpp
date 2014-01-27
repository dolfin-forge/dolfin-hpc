// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLMixedElement.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
MixedElement::MixedElement(
    FiniteElementBaseList const& elements ) :
    FiniteElementBase(ElementList::Mixed, get_cell(elements),
                         get_degree_max(elements)),
    sub_elements_(elements)
{
  // Check mixed finite element definition
  if(elements.size() < 2)
  {
    error("A mixed element should contain more than one subelement");
  }

  // Create string representation
  std::stringstream ssrepr;
  std::stringstream ssstr;
  ssrepr << "MixedElement(*[";
  ssstr << "<Mixed element: (";

  FiniteElementBaseList::const_iterator it = sub_elements_.begin();
  uint value_size_sum = (*it)->value_shape().size();
  ssrepr << (*it)->repr();
  ssstr << (*it)->str();
  for ( ++it ; it != sub_elements_.end(); ++it)
  {
    ssrepr << ", " << (*it)->repr();
    ssstr << ", " << (*it)->str();

  }

  ssrepr << "], **{'value_shape': " << value_shape_.str() << " })";
  ssstr <<  ")>";
  repr_ = ssrepr.str();
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
MixedElement::~MixedElement()
{
}

//-----------------------------------------------------------------------------
bool const MixedElement::is_cellwise_constant() const
{
  bool ret = true;
  for ( FiniteElementBaseList::const_iterator it = sub_elements_.begin();
        it != sub_elements_.end(); ++it )
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<uint, uint> const MixedElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const MixedElement::extract_subelement_component(
    ValueArray const& i) const
{
  check_component(i);

  ValueArray k;
  ValueArray j;
  /*
  if(value_shape_.size() == 1)
  {
    j.insert(i.begin(),i.end());
    for()
    {

    }
  }
  else
  {
    k = i[0];
    j = i.pop
  }

         # Select between indexing modes
         if len(self.value_shape()) == 1:
             # Indexing into a long vector of flattened subelement shapes
             j, = i

             # Find subelement for this index
             for k, e in enumerate(self._sub_elements):
                 sh = e.value_shape()
                 si = product(sh)
                 if j < si:
                     break
                 j -= si
             ufl_assert(j >= 0, "Moved past last value component!")

             # Convert index into a shape tuple
             j = index_to_component(j, sh)
         else:
             # Indexing into a multidimensional tensor
             # where subelement index is first axis
             k = i[0]
             ufl_assert(k < len(self._sub_elements),
                        "Illegal component index (dimension %d)." % k)
             j = i[1:]
  */
  return std::pair<ValueArray, ValueArray>(k,j);
}

//-----------------------------------------------------------------------------
std::pair<uint, FiniteElementBase const * const> const MixedElement::extract_component(ValueArray const& i) const
{
  ValueArray subidx;
  subidx.insert(subidx.begin(),i.begin()+1,i.end());
  return sub_elements_[i[0]]->extract_component(subidx);
}

//-----------------------------------------------------------------------------
uint const MixedElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& MixedElement::sub_elements() const
{
  return sub_elements_;
}

}

