//
//  Version: $Id$
//
#ifndef _VFF_VA_IMPL_H_
#define _VFF_VA_IMPL_H_

#include <opt/debug.h>
#include <opt/tinyxml.h>
 
namespace LaDa
{
  namespace Vff
  {

    template< class T_VFFBASE >
    bool VABase<T_VFFBASE> :: Load( const TiXmlElement &_node )
    {
      namespace bfs = boost::filesystem;
      const TiXmlElement* parent
        = opt::find_node( _node, "Functional", "type", "vff" );
      __DOASSERT( not parent, 
                  "Could not find <Functional type=\"vff\"> in input\n"; )
      bfs::path path;
      TiXmlDocument doc;
      if(  parent->Attribute( "filename" ) )
      {
        path = Print::reformat_home( parent->Attribute( "filename" ) );
        __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
        opt::read_xmlfile( path, doc );
        __DOASSERT( not doc.FirstChild( "Job" ),
                    "Root tag <Job> does not exist in " + path.string() + ".\n" )
        parent = opt::find_node( *doc.FirstChildElement( "Job" ),
                                 "Functional", "type", "vff" );
        __DOASSERT( not parent, 
                    "Could not find <Functional type=\"vff\"> in input\n"; )
      }
      // Load base
      if( not t_VffBase :: Load( _node ) ) return false;

      if( minimizer.Load( *parent ) ) return true;
      return minimizer.Load( *parent->Parent()->ToElement() );
    }


    template< class T_VFFBASE > typename VABase<T_VFFBASE> :: t_Type 
      VABase<T_VFFBASE> :: evaluate()
      {
        // no minimization required if variables is empty.
        typename t_VffBase :: t_Arg arg;
        t_VffBase :: init( arg );
          
        if( arg.size() ) minimizer( *( (t_VffBase*) this), arg );
     
        t_VffBase :: structure.energy = t_VffBase::energy();

        return t_VffBase::structure.energy;
      }

//   template< class T_VFFBASE > typename VABase<T_VFFBASE> :: t_Type 
//     VABase<T_VFFBASE> :: evaluate_one_gradient( types::t_unsigned _pos )
//     {
//       __ASSERT( _pos > va_vars.size(),
//                 "Requesting out-of-range gradient.\n")
//     
//       typename t_Centers :: iterator i_center = centers.begin();
//       typename t_Centers :: iterator i_center_end = centers.end();
//       for(++_pos; _pos and i_center != i_center_end; ++i_center )
//         if( not (i_center->Origin().freeze & t_Atom::FREEZE_T) ) --_pos;
//     
//       t_Type result = functionals[i_center->kind()].evaluate( *i_center );
//       i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//       result -= functionals[i_center->kind()].evaluate( *i_center );
//       result /= t_Type(2);
//       i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//     
//       return result;
//     }
//
//   template< class T_VFFBASE >
//   void VABase<T_VFFBASE> :: evaluate_gradient( t_Type * _grad )
//   {
//     t_Type* i_grad = _grad;
//     typename t_Centers :: iterator i_center = centers.begin();
//     typename t_Centers :: iterator i_center_end = centers.end();
//     for(; i_center != i_center_end; ++i_center )
//     {
//       if( i_center->Origin().freeze & t_Atom::FREEZE_T ) continue;
//
//       *i_grad = functionals[i_center->kind()].evaluate( *i_center );
//       i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//       *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
//       *i_grad /= t_Type(2);
//       i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//
//       ++i_grad;
//     } 
//
//   }
//
//   template< class T_VFFBASE > typename VABase<T_VFFBASE> :: t_Type 
//     VABase<T_VFFBASE> :: evaluate_with_gradient( t_Type * _grad )
//     {
//       t_Type result(0);
//       t_Type* i_grad = _grad;
//       typename t_Centers :: iterator i_center = centers.begin();
//       typename t_Centers :: iterator i_center_end = centers.end();
//       for(; i_center != i_center_end; ++i_center )
//       {
//         if( i_center->Origin().freeze & t_Atom::FREEZE_T ) 
//         {
//           result += functionals[i_center->kind()].evaluate( *i_center );
//           continue;
//         }
//     
//         *i_grad = functionals[i_center->kind()].evaluate( *i_center );
//         result += *i_grad;
//         i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//         *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
//         *i_grad /= t_Type(2);
//         i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
//     
//         ++i_grad;
//       } 
//     
//       return result;
//     }
//
    template< class T_VFFBASE > 
      inline bool VABase<T_VFFBASE> :: init( bool _redocenters )
      {
        if ( _redocenters and ( not t_VffBase::initialize_centers() ) ) return false;
//        if ( _redocenters and ( not t_VffBase::build_tree() ) ) return false;
        return t_VABase::init();
      }

  } // namespace VFF
} // namespace LaDA
#endif // _VFF_VA_IMPL_H_
