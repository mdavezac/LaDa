//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "translation.h"
#include "exceptions.h"


namespace LaDa
{

  namespace enumeration
  {
     t_uint Translation::operator()(t_uint _x, FlavorBase const &_flavorbase) const
     {
       if(card_ < 2) return _x;
#      ifdef LADA_DEBUG
         if(_flavorbase.size() != card_) BOOST_THROW_EXCEPTION( argument_error() );
         if(_x >= _flavorbase.back() * _flavorbase[1]) BOOST_THROW_EXCEPTION( integer_too_large());
#      endif
       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       for(size_t d(0); d < nsites_; ++d)
       {
         atat::iVector3d g(0,0,0);
         for(g(0)=0; g(0) < smith_(0); ++g(0))
         {
           for(g(1)=0; g(1) < smith_(1); ++g(1))
           {
             for(g(2)=0; g(2) < smith_(2); ++g(2), ++i_flavor)
             {
               atat::iVector3d t
               (
                 (g(0) + translation_(0)) % smith_(0),
                 (g(1) + translation_(1)) % smith_(1),
                 (g(2) + translation_(2)) % smith_(2)
               );
               if( t(0) < 0 ) t(0) += smith_(0);
               if( t(1) < 0 ) t(1) += smith_(1);
               if( t(2) < 0 ) t(2) += smith_(2);
               
               t_uint const flavor( _x / (*i_flavor) );
               _x %= (*i_flavor);

               result += flavor * _flavorbase[get_index(d, t, smith_, card_)];
             } // g2
           } // g1
         } // g0
       } // d
       return result;
     }


    boost::shared_ptr< std::vector<Translation> > 
      create_translations(atat::iVector3d const &_smith, size_t _nsites)
      {
        boost::shared_ptr< std::vector<Translation> > result( new std::vector<Translation> );
        atat::iVector3d g(0,0,0);
        for(; g(0) < _smith(0); ++g(0))
        {
          for(g(1)=0; g(1) < _smith(1); ++g(1))
          {
            for(g(2)=0; g(2) < _smith(2); ++g(2))
              if( g(0) != 0 or g(1) != 0 or g(2) != 0 )
                result->push_back( Translation(g, _smith, _nsites) );
          }
        }
        return result;
      }

  } // namespace enumeration.
} // namespace LaDa
