//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_POLICICES_H_
#define _DARWIN_ALLOY_LAYERS_POLICICES_H_

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
          //! \brief Initializes the LayerDepth functor.
          //! \note This member is called by the evaluator. It does not need to
          //!       exist in all translate policies, but if it does, it cannot be
          //!       static.
          template< class T_THIS >
            void init( const T_THIS &_this )
             { Translate<T_OBJECT>::layerdepth.set( _this.structure.cell ); }

          using BitString::Translate< T_OBJECT > :: translate;

          protected:
            //! A layer depth operator.
            static Crystal::LayerDepth layerdepth;
        };
      
      //! \brief A function to easily create random individuals.
      //! \details The object of \a _indiv is translated from a random structure.
      template< class T_INDIVIDUAL, class T_TRANSLATE >
        bool initialize( T_INDIVIDUAL &_indiv,
                         Crystal::Structure& _structure,
                         T_TRANSLATE _translate );
                                   

      template< class T_OBJECT >
        Crystal::LayerDepth Translate<T_OBJECT> :: layerdepth;


      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const Crystal::Structure& _structure,
                                               t_Object& _object )
        {
          // Clears current object.
          typename t_Object :: t_Container& container = _object.Container();
          container.clear();

          typedef typename Crystal::Structure::t_Atoms::const_iterator t_iatom;
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          types::t_real last_depth( Translate<T_OBJECT>::layerdepth( i_atom->pos ) );
          types::t_unsigned layer_size(0);
          types::t_real c(0);
          for( ; i_atom != i_atom_end; )
          {
            types::t_real new_depth;
            do
            {
              c += i_atom->type > 0e0 ? 1e0: -1e0;
              ++layer_size;
              ++i_atom;
              new_depth = Translate<T_OBJECT>::layerdepth( i_atom->pos );
            }
            while( Fuzzy::eq( last_depth, new_depth ) and i_atom != i_atom_end );

            // New layer has been reached
            container.push_back( c / (types::t_real) layer_size );

            // Reinitializing layer discovery.
            layer_size = 0;
            last_depth = new_depth;
            c = 0e0;
          }
        }

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const t_Object& _object,
                                               Crystal::Structure& _structure ) 
        {
          types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& )
              = &eo::random<types::t_unsigned>;
          typedef typename Crystal::Structure::t_Atoms::iterator t_iatom;
          typedef typename t_Object :: t_Container :: const_iterator t_ivar;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          t_iatom i_atom_begin( i_atom );
          t_ivar i_var = _object.Container().begin();
          __DODEBUGCODE( t_ivar i_var_end = _object.Container().end(); )
          types::t_real last_depth( Translate<T_OBJECT>::layerdepth( i_atom->pos ) );
          types::t_unsigned layer_size(0);
          for( ; i_atom != i_atom_end __DODEBUGCODE( and i_var != i_var_end ); )
          {
            types::t_real new_depth;
            do
            {
              ++layer_size;
              ++i_atom;
              new_depth = Translate<T_OBJECT>::layerdepth( i_atom->pos );
            }
            while( Fuzzy::eq( last_depth, new_depth ) and i_atom != i_atom_end );

            // New layer has been reached
            // Creates an array with approximate concentration, shuffles, then
            // assigns to atoms.
            const bool is_positive( *i_var > 0 );
            const bool is_negative(not is_positive);
            const types::t_real nb_layers_real
            (
              is_positive ? ( 1e0 - (*i_var) ) * 0.5e0 * types::t_real( layer_size ):
                            ( 1e0 + (*i_var) ) * 0.5e0 * types::t_real( layer_size )
            );
            const size_t nb_layers_size_t( nb_layers_real );
            std::vector<bool> arrangment(layer_size,  is_positive );
            std::fill
            ( 
              arrangment.begin(),
              arrangment.begin() + nb_layers_size_t,
              is_negative
            );
            std::random_shuffle( arrangment.begin(), arrangment.end(), ptr_to_rng );
            std::vector<bool> :: const_iterator i_bool = arrangment.end(); 
            for(; i_atom_begin != i_atom; ++i_atom_begin, ++i_bool )
              i_atom_begin->type = *i_bool ? 1e0: -1e0;

            // Reinitializing layer discovery.
            layer_size = 0;
            last_depth = new_depth;
            i_atom_begin = i_atom;
            ++i_var;
          }
          __ASSERT( i_atom != i_atom_end, "Structure smaller than object.\n" )
          __ASSERT( i_var != i_var_end, "Object smaller than structure.\n" )
        }
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa

#endif
