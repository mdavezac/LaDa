//
//  Version: $Id$
//
#ifndef _LADA_MINIMIZER_ANY_H_
#define _LADA_MINIMIZER_ANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/mpl/vector.hpp>

#include <opt/tinyxml.h>

namespace LaDa
{
  namespace Minimizer
  {
    //! Wraps around other minimimizers to provided single Load/Launch.
    template< class T_MPL_VECTOR >
      class Variant
      {
        public:
          //! list of minimizers.
          typedef T_MPL_VECTOR t_Minimizers;

          //! Constructor.
          Variant() {}
          //! Copy Constructor.
          Variant( const Variant& _c ) : minimizer( _c.minimizer ) {}
  
  
          //! Calls minimizer.
          template< class T_FUNCTOR >
            typename T_FUNCTOR :: t_Return 
              operator()( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b ) const;
          
          //! Load Frpr or Gsl minimizer from XML.
          bool Load( const TiXmlElement& _node );
  
        protected:
          //! A points to a Gsl minimizer.
          typename boost::make_variant_over< t_Minimizers > :: type minimizer;
      };

    namespace details
    {
      template< class T_FUNCTION >
        struct OpVisitor : public boost::static_visitor< typename T_FUNCTION :: t_Return > 
        {
          const T_FUNCTION &func_;
          typename T_FUNCTION :: t_Arg &arg_;
        
          OpVisitor   ( const T_FUNCTION &_func, typename T_FUNCTION :: t_Arg &_arg )
                    : func_( _func ), arg_( _arg ) {}
          OpVisitor( const OpVisitor &_c ) : func_( _c.func_ ), arg_( _c.arg_ ) {}

          template< class T_MINIMIZER >
            typename T_FUNCTION :: t_Return operator()( const T_MINIMIZER& _minimizer ) const
              { return _minimizer( func_, arg_ ); }
        };


      //! Registers a minimizer for variant loading.
      template< class T_MINIMIZER > struct MinimizerVariantType;

      template< class T_FIRST, class T_LAST, class T_VARIANT >
        struct LoadVariant
        {
          protected:
            struct load 
            {
              template< class T_FUNCTOR >
                static bool apply( T_FUNCTOR& _func, const TiXmlElement &_node )
                {
                  typedef typename boost::mpl::deref<T_FIRST> :: type t_Type;
                  if( not _func.Load( _node ) ) return false;
                  std::cerr <<   "Successfuly loaded "
                               + MinimizerVariantType<t_Type>::name
                               + " minimizer.\n";
                  return true;
                }
            };
      
          public:
            bool static apply( T_VARIANT &_variant, const TiXmlElement& _node );
        };
      template< class T_FIRST, class T_LAST, class T_VARIANT >
        bool LoadVariant<T_FIRST, T_LAST, T_VARIANT> :: 
          apply( T_VARIANT &_variant, const TiXmlElement& _node )
          {
            typedef typename boost::mpl::deref<T_FIRST> :: type t_Type;
            t_Type minimizer;
            if( not load::apply( minimizer, _node ) )
              return LoadVariant
                     < 
                       typename boost::mpl::next< T_FIRST > :: type, 
                       T_LAST, 
                       T_VARIANT 
                     > :: apply( _variant, _node );
            _variant = minimizer;
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
          return LoadVariant
                 < 
                   typename boost::mpl::begin< T_SEQUENCE > :: type, 
                   typename boost::mpl::end< T_SEQUENCE > :: type, 
                   typename boost::make_variant_over< T_SEQUENCE > :: type 
                 > :: apply( _variant, _node );
        }

    }

    template< class T_MINIMIZERS > template< class T_FUNCTOR >
      typename T_FUNCTOR :: t_Return 
        Variant<T_MINIMIZERS> :: operator()( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b ) const
        {
          typedef typename T_FUNCTOR :: t_Return t_Return;
          typedef typename T_FUNCTOR :: t_Arg t_Arg;
          return boost::apply_visitor( details::OpVisitor<T_FUNCTOR>( _a, _b ), minimizer );
        }


  template< class T_MINIMIZERS > 
    bool Variant<T_MINIMIZERS>::Load( const TiXmlElement& _node )
    {
      const TiXmlElement *node = opt::find_node( _node, "Minimizer" );
      for(; node; node = node->NextSiblingElement( "Minimizer" ) )
        if( details::load_variant<t_Minimizers>( minimizer, *node ) ) return true;

      return false;
    }
  }
# define LADA_REGISTER_MINIMIZER_VARIANT_HEADER( a, b ) \
    namespace details \
    {  \
      template<> \
        struct MinimizerVariantType<a> \
        { \
          typedef a type; \
          static const std::string name; \
        }; \
    }
# define LADA_REGISTER_MINIMIZER_VARIANT_SOURCE( a, b ) \
    namespace details \
    {  \
      const std::string MinimizerVariantType<a> :: name = b; \
    }
} 
#endif
