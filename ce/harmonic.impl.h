//
//  Version: $Id$
//

namespace LaDa
{
  namespace CE
  {
    namespace ConstituentStrain
    {
      namespace Harmonic
      {

        template< class T_DERIVED > inline 
          types::t_real Base<T_DERIVED> :: evaluate(const atat::rVector3d &_k) const
          {
            return (  exp( -norm2(_k) * attenuation ) * 
                   (*static_cast<const t_Derived*>(this))( _k ) ); 
          }

        template< class T_DERIVED > inline 
        types::t_real Base<T_DERIVED> :: evaluate(const types::t_real _x,
                                                  const atat::rVector3d &_k) const
        {
          return (   interpolation.evaluate(_x) 
                   * exp( -norm2(_k) * attenuation )
                   * (*static_cast<const t_Derived*>(this))( _k ) );
        }
        template< class T_DERIVED > inline 
        types::t_real Base<T_DERIVED> :: evaluate_gradient(const types::t_real _x,
                                                           const atat::rVector3d &_k) const
        {
          return (   interpolation.evaluate_gradient(_x) 
                   * exp( -norm2(_k) *  attenuation )
                   * (*static_cast<const t_Derived*>(this))( _k ) );
        }
        template< class T_DERIVED > inline 
        types::t_real Base<T_DERIVED> :: evaluate_with_gradient(const types::t_real _x, 
                                                                const atat::rVector3d &_k,
                                                                types::t_real &_grad) const
        {
          types::t_real factor = exp( -norm2(_k) * attenuation )
                          * (*static_cast<const t_Derived*>(this))( _k );
          types::t_real result =   interpolation.evaluate_with_gradient(_x, _grad) 
                          * factor ;
          _grad *= factor;
          return result;
        }

        template< class T_DERIVED > 
        bool Base<T_DERIVED> :: Load (const TiXmlElement &_element) 
        {
          const TiXmlElement *child;
          int i=1;

          if( _element.Attribute("type" ) )
          {
            std::string str = _element.Attribute("type");
            if( str.compare( t_Derived::type ) ) return false;
          }
                  
          __DOASSERT( not _element.Attribute("rank", &i),
                      "Harmonic has no rank on input.\n" )
          __DOASSERT( i < 0 or i > t_Derived::maxrank,
                      "rank of harmonic is out of range on input.\n" )

          rank = types::t_unsigned(abs(i));
          child = _element.FirstChildElement( "Point" );
          clear(); // clear interpolation
          __DOASSERT( not child, "No interpolation points found for harmonic on input.\n" )
          for ( ; child; child = child->NextSiblingElement( "Point" ) )
            __DOASSERT( not interpolation.Load( *child ),
                        "Error while loading harmonic interpolation point.\n" )

          return true;
        }


      } // namespace Harmonic
    } // namespace ConstituentStrain
  } // namespace CE
    
} // namespace LaDa
