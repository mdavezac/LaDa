//
//  Version: $Id$
//
#include <iostream>
#include <sstream>
#include <string>
#include <functional>
#include <cmath>
#include <math.h>
#include <string>
#include <stdio.h>


#include "analysis/analyze_code.h"
#include "analysis/call_log.h"

#include "polynome.h"

namespace VA_CE {

  types::t_unsigned Polynome::nb_eval           = 0;
  types::t_unsigned Polynome::nb_eval_grad      = 0;
  types::t_unsigned Polynome::nb_eval_with_grad = 0;


  // resizes variables vector such that it matches exactly the terms
  // in the monomes vector
  void Polynome :: resize_variables()
  {
    types::t_int m1 = 0, m2 = 0;
    std::list< function::Monome<> > :: iterator i_monome = monomes.begin();
    std::list< function::Monome<> > :: iterator i_last   = monomes.end();
    for( ; i_monome != i_last ; i_monome++)
    {
      m1 = i_monome->max_term();
      if ( m1 > m2 )
        m2 = m1;
    }
    m2++;
    if ( m2 > types::t_int(variables->size()) )
      variables->resize(m2, 0);
    else if ( m2 < types::t_int(variables->size()) )
      variables->resize(m2);
  }


  void Polynome :: print_xml( TiXmlElement* const node ) const
  {
    TiXmlElement *parent;

    parent = new TiXmlElement( "Polynomial" );
    node->LinkEndChild( parent );

    if ( monomes.size() == 0 )
      return; 

    std::list< function::Monome<> > :: const_iterator i_monome = monomes.begin(); 
    std::list< function::Monome<> > :: const_iterator i_end = monomes.end(); 
    for( ;  i_monome != i_end; ++i_monome )
      print_xml(*i_monome, parent);
  }

  bool Polynome :: Load( TiXmlElement* const node )
  {
    TiXmlElement *child;
    function::Monome<> monome;
    monomes.clear();

    if ( not node ) 
      return false;
    child = node->FirstChildElement("Polynomial");
    if ( not child )
      return false;

    child = child->FirstChildElement("Monomial");
    for( ;  child; child  =  child->NextSiblingElement("Monomial") )
    {
      if ( not Load(monome, child) )
        return false;
      monomes.push_back(monome);
    }

    if ( not variables )
      variables = new t_Container;
      
    return true;
  }


  

  #ifdef _DEBUG_LADA_ 
    bool Polynome :: is_consistent() const
    {
      return true;
//     if ( not variables )
//     {
//       std::cerr << " variables were not initialized in Polynome object"
//                 << std::endl;
//       exit(0);
//     }
//
//     bool is_valid = true;
//     if ( variables->size() == 0 || monomes.size() == 0)
//     {
//       std::cerr << " Variables or Monomes are not initialized "
//                 << std::endl;
//       is_valid = false;
//     }
//      
//     std::list< function::Monome<> > :: const_iterator i_monome = monomes.begin();
//     std::list< function::Monome<> > :: const_iterator i_monome_last = monomes.end();
//     for ( ; i_monome != monomes.end(); i_monome++ ) // monome loop
//     {
//       std::list<types::t_int> :: const_iterator i_term = i_monome->terms.begin();
//       std::list<types::t_int> :: const_iterator i_term_last = i_monome->terms.end();
//     
//       for( ; i_term != i_term_last; i_term++)
//         if ( variables->size() <= abs(*i_term) )
//         {
//           std::cerr << " term of monome is larger than the number "
//                     << "of variables in polynome "
//                     << std::endl;
//           is_valid = false;
//         }
//     
//     } // end of loop over monomes
//
//     return is_valid;
    }
  #endif // _DEBUG_LADA_ 
       
  #ifdef _DEBUG_LADA_FAKE_POLY_
    void Polynome::fake_poly()
    {
      variables->clear();
      variables->push_back(0.0);
      variables->push_back(0.0);
      monomes.clear();
      function::Monome<> monome;
      monome.coefficient = 1;
      monome.add_term(0);
      monome.add_term(1);
      add(monome);
    }
  #endif // _DEBUG_LADA_FAKE_POLY_

       
  void Polynome :: print_xml( const function::Monome<> &_m, TiXmlElement* const node ) const
  {
    TiXmlElement *parent;

    parent = new TiXmlElement( "Monomial" );
    node->LinkEndChild( parent );
    parent->SetDoubleAttribute( "coefficient", _m.get_coefficient() );

    std::ostringstream stream; 
    _m.print_terms(stream);
    parent->SetAttribute( "terms", stream.str().c_str() );
  }

  bool Polynome :: Load( function::Monome<> &_m, TiXmlElement* const node )
  {
//   _m.terms.clear();
//
//   if ( not node ) 
//     return false;
//
//   _m.coefficient = 1.0;
//   if ( not node->Attribute("coefficient", &coefficient) )
//     return false;
//
//   const char* cchar = node->Attribute("terms");
//   if ( not cchar ) // zero order term;
//     return true; 
//
//   // order >0 terms
//   std::istringstream istr ( cchar );
//   
//   if ( istr.eof() )
//     return false;
//
//   types::t_int d;
//   while ( not istr.eof() )
//   { 
//     istr >> d;
//     _m.terms.push_back(d);
//   }
//
    return true;
  }

} // namespace VA_CE
