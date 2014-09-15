// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-03
// Last changed: 2014-02-03

#include <dolfin/fem/FiniteElementSpace.h>

#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofMap.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(
    Mesh& mesh, std::string const& finite_element_signature,
    std::string const& dof_map_signature) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(finite_element_signature),
    dof_map_(DofMapCache::instance().acquire_dofmap(mesh, dof_map_signature))
#if ENABLE_UFL
                                                    ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(finite_element_signature)))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, std::string const& signature) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(signature),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(signature)))
#if ENABLE_UFL
                                           ,
    ufl_(ufl::FiniteElementBase::create(ufl::Object::repr_t(signature)))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, Form& form, uint const i) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(mesh.type(), form, i),
    dof_map_(DofMapCache::instance().acquire_dofmap(mesh, form, i))
#if ENABLE_UFL
                                                    ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh,
                                       ufc::finite_element& finite_element,
                                       bool const finite_element_local) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(finite_element, finite_element_local),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(finite_element_.signature())))
#if ENABLE_UFL
                                           ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
#endif
{
}

#if ENABLE_UFL
//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(
    Mesh& mesh, ufl::FiniteElementBase const& finite_element) :
    mesh_(mesh),
    cell_(mesh, 0),
    finite_element_(finite_element),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(finite_element_.signature()))),
    ufl_(&finite_element)
{
}
#endif

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& other) :
    mesh_(other.mesh()),
    cell_(other.cell()),
    finite_element_(other.element()),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(other.mesh(),
                                               other.dofmap().signature()))
#if ENABLE_UFL
                                               ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& space,
                                       uint const& i) :
    mesh_(space.mesh()),
    cell_(space.cell()),
    finite_element_(*space.element().create_sub_element(i), true),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            space.mesh(),
            DofMap::dofmap_signature(finite_element_.signature())))
#if ENABLE_UFL
                                     ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::~FiniteElementSpace()
{
  DofMapCache::instance().release_dofmap(dof_map_);
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::operator ==(FiniteElementSpace const& other) const
{
  return (this->mesh() == other.mesh())
      && (std::strcmp(this->element().signature(), other.element().signature())
          == 0)
      && (std::strcmp(this->dofmap().signature(), other.dofmap().signature())
          == 0);
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
  return (mesh_.global_numCells() * dof_map_.local_dimension())
      == dof_map_.global_dimension();
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_cellwise_constant() const
{
  return is_cellwise_defined()
      && (dof_map_.local_dimension() == finite_element_.value_dimension(0));
}

//-----------------------------------------------------------------------------

}
/* namespace icorne */
