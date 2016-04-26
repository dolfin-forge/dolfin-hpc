// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLMixedElement.h>

#include <dolfin/ufl/UFLarray.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
MixedElement::MixedElement(FiniteElementSpace::List const& elements) :
    FiniteElementSpace("MixedElement"),
    sub_elements_(cloneSubElementsList(elements)),
    family_(Family::Mixed),
    cell_(get_cell(elements)),
    degree_(get_degree_max(elements)),
    value_shape_(1, get_value_size(elements))
{
  createReprStr();
}

//-----------------------------------------------------------------------------
MixedElement::MixedElement(repr_t const& repr) :
    FiniteElementSpace("MixedElement", repr),
    sub_elements_(createSubElementsList(repr)),
    family_(Family::Mixed),
    cell_(get_cell(sub_elements_)),
    degree_(get_degree_max(sub_elements_)),
    value_shape_(1, get_value_size(sub_elements_))
{
  createReprStr();
}

//-----------------------------------------------------------------------------
MixedElement::~MixedElement()
{
  // Cleanup sub elements which were either cloned or created from the
  // representation string
  for (FiniteElementSpace::List::const_iterator it = sub_elements_.begin();
      it != sub_elements_.end(); ++it)
  {
    delete (*it);
  }
}

//-----------------------------------------------------------------------------
Family const& MixedElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
Family::Type MixedElement::metatype() const
{
  return Family::Mixed;
}

//-----------------------------------------------------------------------------
Cell const& MixedElement::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& MixedElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueArray const& MixedElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool MixedElement::is_cellwise_constant() const
{
  bool ret = true;
  for (List::const_iterator it = sub_elements_.begin();
      it != sub_elements_.end(); ++it)
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const& MixedElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> MixedElement::extract_subelement_component(
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
  return std::pair<ValueArray, ValueArray>(k, j);
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementSpace const *> MixedElement::extract_component(
    ValueArray const& i) const
{
  ValueArray subidx;
  subidx.insert(subidx.begin(), i.begin() + 1, i.end());
  return sub_elements_[i[0]]->extract_component(subidx);
}

//-----------------------------------------------------------------------------
dolfin::uint MixedElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementSpace::List const& MixedElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& MixedElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& MixedElement::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
FiniteElementSpace::List MixedElement::cloneSubElementsList(
    FiniteElementSpace::List const& elements)
{
  FiniteElementSpace::List subelms;
  //FIXME: Improve performance by implementing a cloning method.
  for (FiniteElementSpace::List::const_iterator it = elements.begin();
      it != elements.end(); ++it)
  {
    subelms.push_back(FiniteElementSpace::create((*it)->repr()));
  }
  dolfin_assert(subelms.size() > 1);
  return subelms;
}

//-----------------------------------------------------------------------------
FiniteElementSpace::List MixedElement::createSubElementsList(repr_t const& repr)
{
  List subelms;
  // Unpack array then split arguments
  std::vector<Object::repr_t> const subreprs = make_args_repr(
      array::unpack(arg(0)), true);
  for (std::vector<Object::repr_t>::const_iterator it = subreprs.begin();
      it != subreprs.end(); ++it)
  {
    subelms.push_back(FiniteElementSpace::create(*it));
  }
  dolfin_assert(subelms.size() > 1);
  return subelms;
}

//-----------------------------------------------------------------------------
void MixedElement::createReprStr()
{
  // Check mixed finite element definition
  if (sub_elements_.size() < 2)
  {
    error("A mixed element should contain more than one subelement");
  }
  
  // Create string representation
  std::stringstream ssrepr;
  std::stringstream ssstr;
  ssrepr << "MixedElement(*[";
  ssstr << "<Mixed element: (";
  
  List::const_iterator it = sub_elements_.begin();
  ssrepr << (*it)->repr();
  ssstr << (*it)->str();
  for (++it; it != sub_elements_.end(); ++it)
  {
    ssrepr << ", " << (*it)->repr();
    ssstr << ", " << (*it)->str();
    
  }
  
  ssrepr << "], **{'value_shape': " << value_shape_.str() << " })";
  ssstr << ")>";
  repr_ = ssrepr.str();
  str_ = ssstr.str();
}

}

