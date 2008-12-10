//
//  Version: $Id$
//
#ifndef _LADA_GA_PURELAYERS_POLICICES_H_
#define _LADA_GA_PURELAYERS_POLICICES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/signal.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

#include <tinyxml/tinyxml.h>

#include <crystal/layerdepth.h>

#include "../bitstring/translate.h"


namespace LaDa
{
  namespace Crystal { class Structure; }

  namespace GA
  {
    namespace AlloyLayers
    {
      //! \brief Policy to translate an object containing a vector of reals to an
      //!        epitaxial alloy structure.
      //! \details It is expected that the atoms in the Crystal::Structure
      //!          instances are ordered in direction of growth ( \see
      //!          bool Crystal::create_epitaxial_structure() ). 
      template< class T_OBJECT >
        struct Translate : public BitString::Translate< T_OBJECT >
        {
          //! Type of the Object.
          typedef T_OBJECT t_Object;
          //! From objects to Crystal::Structure. 
          static void translate( const t_Object&, Crystal :: Structure& );
          //! From Crystal::Structure to objects.
          static void translate( const Crystal :: Structure&, t_Object& );

          using BitString::Translate< T_OBJECT > :: translate;

          protected:
            //! A layer depth operator.
            static Crystal::LayerDepth layerdepth;
        };
      

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const Crystal::Structure& _structure,
                                               t_Object& _object )
        {
          typename t_Object :: t_Container& container = _object.Container();
          container.clear();

          typedef typename Crystal::Structure::t_Atoms::const_iterator t_iatom;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          for( ; i_atom != i_atom_end; ++i_atom )
            if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T ) 
              container.push_back( i_atom->type > 0 ? 1: 0 );
        }

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const t_Object& _object,
                                               Crystal::Structure& _structure ) 
        {
          typedef typename T_OBJECT :: t_Container :: const_iterator t_ivar;
          typedef typename Crystal::Structure::t_Atoms::iterator t_iatom;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          t_ivar i_var = _object.begin(); 
          for( ; i_atom != i_atom_end; ++i_atom, ++i_var )
            if ( not (i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T) )
              i_atom->type = *i_var > 0 ? 1.e0: -1e0;
        }
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa

#endif
