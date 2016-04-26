// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-03
// Last changed: 2014-02-03

#include <dolfin/fem/FiniteElementSpace.h>

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/SubSystem.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Form& form, uint const i) :
    mesh_(form.dofmaps()[i].mesh()),
    cell_(mesh_, 0),
    finite_element_(mesh_.type(), form, i),
    dof_map_(DofMap::acquire(mesh_, form, i)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, Form& form, uint const i) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(mesh.type(), form, i),
    dof_map_(DofMap::acquire(mesh, form, i)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, ufc::finite_element& element,
                                       ufc::dofmap& dofmap, bool owner) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(element, owner),
    dof_map_(DofMap::acquire(mesh, dofmap, owner)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(finite_element_.signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh,
                                       ufl::FiniteElementSpace const& element) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(element),
    dof_map_(
        DofMap::acquire(
            mesh,
            *ElementLibrary::create_dof_map(
                DofMap::make_signature(finite_element_.signature())),
            true)),
    ufl_(ufl::FiniteElementSpace::create(element.repr()))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& space,
                                       uint const i) :
    mesh_(space.mesh()),
    cell_(space.cell()),
    finite_element_(space.element(), i),
    dof_map_(
        DofMap::acquire(space.mesh(), *space.dofmap().create_sub_dofmap(i),
                        true)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& space,
                                       SubSystem const& sub) :
    mesh_(space.mesh()),
    cell_(space.cell()),
    finite_element_(space.element(), sub),
    dof_map_(
        DofMap::acquire(space.mesh(), *space.dofmap().create_sub_dofmap(sub),
                        true)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& other_mesh,
                                       FiniteElementSpace const& space) :
    mesh_(other_mesh),
    cell_(space.cell()),
    finite_element_(space.element()),
    dof_map_(DofMap::acquire(other_mesh, *space.dofmap().create(), true)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))

{
  if (other_mesh.type().cellType() != space.cell().type())
  {
    error("Provided mesh is incompatible with given space");
  }
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& other) :
    mesh_(other.mesh()),
    cell_(other.cell()),
    finite_element_(other.element()),
    dof_map_(DofMap::acquire(other.mesh(), *other.dofmap().create(), true)),
    ufl_(
        ufl::FiniteElementSpace::create(
            ufl::Object::repr_t(element().signature())))
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::~FiniteElementSpace()
{
  DofMap::release(dof_map_);
  delete ufl_;
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::operator ==(FiniteElementSpace const& other) const
{
  return (this->mesh() == other.mesh()) && (this->element() == other.element())
      && (this->dofmap() == other.dofmap());
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::operator !=(FiniteElementSpace const& other) const
{
  return !(*this == other);
}

//-----------------------------------------------------------------------------
Cell& FiniteElementSpace::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
FiniteElement const& FiniteElementSpace::element() const
{
  return finite_element_;
}

//-----------------------------------------------------------------------------
DofMap const& FiniteElementSpace::dofmap() const
{
  return dof_map_;
}

//-----------------------------------------------------------------------------
Array<FiniteElementSpace *> FiniteElementSpace::flatten() const
{
  Array<FiniteElementSpace *> flt;
  Array<ufc::finite_element const*> flt_fe = finite_element_.flatten();
  Array<ufc::dofmap const*> flt_dm = dof_map_.flatten();

  dolfin_assert(flt_fe.size() == flt_dm.size());
  for (uint s = 0; s < flt_fe.size(); ++s)
  {
    flt.push_back(
        new FiniteElementSpace(mesh_, *flt_fe[s]->create(),
                               *flt_dm[s]->create(), true));
  }
  return flt;
}

//-----------------------------------------------------------------------------
void FiniteElementSpace::disp() const
{
  cout << "FiniteElementSpace" << endl;
  cout << "------------------" << endl;

  // Begin indentation
  begin("");
  cout << "Finite element        : " << this->element().signature() << endl;
  cout << "Dof map               : " << this->dofmap().signature() << endl;
  // End indentation
  end();
  skip();
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_cellwise_defined() const
{
  return (mesh_.num_global_cells() * dof_map_.local_dimension())
      == dof_map_.global_dimension();
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_cellwise_constant() const
{
  return is_cellwise_defined()
      && (dof_map_.local_dimension() == finite_element_.value_size());
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_vertex_based() const
{
  //FIXME: Only a particular case.
  return (this->family() == ufl::Family::CG) && (this->degree() == 1);
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_flattenable() const
{
  //FIXME: Only a particular case.
  uint value_size = 1;
  for(uint i = 0; i < this->element().value_rank(); ++i)
  {
    value_size *= this->element().value_dimension(i);
  }
  return (this->element().flatten().size() == value_size);
}

//-----------------------------------------------------------------------------

}
/* namespace icorne */
