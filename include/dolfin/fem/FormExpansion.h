// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed:

#ifndef __FORM_EXPANSION_H
#define __FORM_EXPANSION_H

#include <dolfin/quadrature/QuadratureRule.h>
#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLElementList.h>
#include <dolfin/ufl/UFLForm.h>


namespace dolfin
{
  class FormExpansion
  {
    public: 

      /// Constructor
      FormExpansion(const ufc::form& form, const QuadratureRule& q);
      
      /// Destructor
      virtual ~FormExpansion();

    private:
      const ufc::form* form;
      const QuadratureRule* quadrature;
      void expand(const std::string& s);
      void parse_form();
      std::vector<std::string> finite_element_signatures;
      std::vector<std::string> expanded_form;
      std::vector<real> values_at_quadrature_points;
  };

  inline FormExpansion::FormExpansion(const ufc::form& f, const QuadratureRule& q)
    : 
      form(&f), quadrature(&q)
  {
    std::cout << "signature    " << form->signature() << std::endl;
    std::string signature(form->signature());
//    std::cout << "signature2   " << signature << std::endl;
//    std::cout << signature.find("Form") << std::endl;
//    std::cout << signature.find("Ich") << std::endl;
//    std::cout << signature.find("MultiIndex") << std::endl;
//    expand(signature);
    ufl::Form ufl_form (signature);
//    parse_form();
  }

  inline FormExpansion::~FormExpansion()
  {}

  inline void FormExpansion::expand(const std::string& s)
  {
    if (s.length() == 0)
    {
      return;
    }
    else
    {
      std::string::const_iterator it;
      std::string work_string;
      for(it=s.begin(); it!=s.end(); ++it)
      {
        work_string += *it;
        if(work_string == "FiniteElement")
        {
//          std::cout << "FE" << std::endl;
          for(std::string::const_iterator jt = ++it; jt!=s.end(); ++jt)
          {
            work_string += *jt;
//            std::cout << "work string:  " << work_string << std::endl;
            if(work_string.find(" None)") != -1)
            {
//              std::cout << "FE push:  " << work_string << std::endl;
              finite_element_signatures.push_back(work_string);
              it = jt;
              break;
            }
          }
//          std::cout << "Form push:  " << work_string << std::endl;
          expanded_form.push_back(work_string);
          work_string.clear();
        }
        else if(work_string == "MultiIndex")
        {
//          std::cout << "MI" << std::endl;
          for(std::string::const_iterator jt = ++it; jt!=s.end(); ++jt)
          {
            work_string += *jt;
//            std::cout << "work string:  " << work_string << std::endl;
            if(work_string.find("))") != -1)
            {
              it = jt;
              break;
            }
          }
//          std::cout << "push:  " << work_string << std::endl;
          expanded_form.push_back(work_string);
          work_string.clear();
        }
        if(work_string == "Measure")
        {
//          std::cout << "M" << std::endl;
          for(std::string::const_iterator jt = ++it; jt!=s.end(); ++jt)
          {
            work_string += *jt;
//            std::cout << "work string:  " << work_string << std::endl;
            if(work_string.find(" None)") != -1)
            {
              it = jt;
              break;
            }
          }
//          std::cout << "M push:  " << work_string << std::endl;
          expanded_form.push_back(work_string);
          work_string.clear();
        }
        else if(*it == '(' || work_string.find(')') != -1)// || *it == ')')    
        {
//          std::cout << "push:  " << work_string << std::endl;
          expanded_form.push_back(work_string);
          work_string.clear();
        }
        if(work_string == ", ")    
        {
//          std::cout << "push:  " << work_string << std::endl;
          expanded_form.push_back(work_string);
          work_string.clear();
        }
      }
      for(unsigned int i=0; i<expanded_form.size(); ++i)
        std::cout << expanded_form[i] << std::endl;
      return;
    }
    return;
  }

  inline void FormExpansion::parse_form()
  {
    std::cout << "parse" << std::endl;
//    values_at_quadrature_points.clear();
//    const std::vector<real*>& q_points = quadrature->get_points();
//    const uint n_q_points = q_points.size();
//    values_at_quadrature_points.resize(n_q_points);

    std::vector<std::string>::iterator it;
    std::vector<std::string>::iterator end_it;
    uint index = 0;
    const uint size = expanded_form.size();
    for(it = expanded_form.begin(); it!=end_it; ++it, ++index)
    {
      std::cout << *it << std::endl; 
      std::cout << expanded_form[size-index-1] << std::endl; 
      if(index == 0 && it->find("Form") == -1)
      {
        dolfin_assert("A form is expected.");
      }
      else
      {
        if(it->find(",")!=-1)
        {
          std::vector<std::string>::iterator node = it;
          node--;
          std::vector<std::string>::iterator leave1 = it;
          std::vector<std::string>::iterator leave2 = it;
          leave2++, leave2++;

          std::cout << "node contains " << *node << std::endl;
          std::cout << "leave1 contains " << *leave1 << std::endl;
          std::cout << "leave2 contains " << *leave2 << std::endl;

          if(node->find("Argument")!=-1)
          {
            if(leave1->find("FiniteElement")==-1)
              dolfin_assert("No FiniteElement given in Argument.");
            else
            {
//              ufl::ElementList::FamilyType family;
              ufl::Cell* cell;
              uint degree;
//              parse_finite_element(*leave1, family, cell, degree);
            }
          }
        }
      }
    }
  }

//  inline void FormExpansion::parse_finite_element(std::string const& s, 
//      ufl::ElementList::FamilyType& family, ufl::Cell* cell, uint & degree) const
//  {
//    std::string::const_iterator it;
//    std::string work_string;
//    for(it=s.begin(); it!=s.end(); ++it)
//    {
//      work_string += *it;
//      if(work_string == "FiniteElement('")
//        work_string.clear();
//      
//      if(work_string == "Lagrange', " || work_string == "CG', ")
//      {
//        family = ufl::ElementList::CG;
//        work_string.clear();
//      }
//
//      if(work_string == "Cell('")
//        work_string.clear();
//
//      ufl::Domain domain;//(ufl::Domain::Type::triangle);
//      if(work_string == "triangle', ")
//      {
//        domain(ufl::Domain::Type::triangle);
//        work_string.clear();
//      }
//
//
//      if(work_string == "Space(")
//      {
//        std::string::const_iterator jt = ++it;
//        std::cout << jt << std::endl;
//        work_string.clear();      
//      }
//    }
//}
}
#endif
