//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "any.h"

namespace LaDa
{
  namespace Minimizer
  {
    template< class T_FIRST, class T_LAST, class T_VARIANT >
      struct LoadVariant
      {
        protected:
          struct load 
          {
            template< class T_FUNCTOR >
              static bool apply( T_FUNCTOR& _func, const TiXmlElement &_node )
              {
                if( not _func.Load( _node ) ) return false;
                std::cerr << "Successfuly loaded " + name(_func) + " minimizer.\n";
                return true;
              }
            const std::string name( Frpr& _func ) { return "original VFF"; }
            const std::string name( Gsl& _func ) { return "GSL"; }
            const std::string name( Decoupled& _func ) { return "decoupled"; }
          };

        public:
          bool static apply( T_VARIANT &_variant, const TiXmlElement& _node );
      };
    template< class T_FIRST, class T_LAST, class T_VARIANT >
      bool LoadVariant<T_FIRST, T_LAST, T_VARIANT> :: 
        apply( T_VARIANT &_variant, const TiXmlElement& _node )
        {
//         typedef typename boost::mpl::deref<T_FIRST> :: type t_Type;
//        t_Type minimizer;
//         if( not load::apply( minimizer, _node ) )
//           LoadVariant
//           < 
//             typename boost::mpl::next< T_FIRST > :: type, 
//             T_LAST, 
//             T_VARIANT 
//           > :: apply( _variant, _node );
//         _variant = minimizer;
          return true;
        }

    template< class T_LAST, class T_VARIANT >
      struct LoadVariant< T_LAST, T_LAST, T_VARIANT >
      {
        static bool apply( T_VARIANT&, const TiXmlElement &);
      };
    template< class T_LAST, class T_VARIANT >
      bool LoadVariant<T_LAST, T_LAST, T_VARIANT> :: 
        apply( T_VARIANT &_variant, const TiXmlElement& _node ) { return false; }

    template< class T_SEQUENCE >
      bool load_variant( typename boost::make_variant_over< T_SEQUENCE > :: type &_variant, 
                         const TiXmlElement &_node )
      {
        return true;
//       return LoadVariant
//              < 
//                boost::mpl::first< T_SEQUENCE >, 
//                boost::mpl::end< T_SEQUENCE >, 
//                typename boost::make_variant_over< T_SEQUENCE > :: type 
//              > :: apply( _variant, _node );
      }


    bool Any::Load( const TiXmlElement& _node )
    {
      const TiXmlElement *node = opt::find_node( _node, "Minimizer" );
      for(; node; node = node->NextSiblingElement( "Minimizer" ) )
        if( load_variant<t_Minimizers>( minimizer, *node ) ) return true;

      return false;
    }

  }
} 
