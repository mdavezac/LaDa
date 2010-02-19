#ifndef LADA_OPT_TUPLES_SERIALIZE_H
#define LADA_OPT_TUPLES_SERIALIZE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/serialization/nvp.hpp>
// taken from http://uint32t.blogspot.com/2008/03/serializing-boosttuple-using.html
namespace boost { namespace serialization {

#   define GENERATE_TUPLE_SERIALIZE(z,nargs,unused)                            \
      template< typename Archive, BOOST_PP_ENUM_PARAMS(nargs,typename T) > \
      void serialize(Archive & ar,                                        \
                     boost::tuple< BOOST_PP_ENUM_PARAMS(nargs,T) > & t,   \
                     const unsigned int version)                          \
      {                                                                   \
        ar & boost::serialization::make_nvp("head",t.head);               \
        ar & boost::serialization::make_nvp("tail",t.tail);               \
      }


      BOOST_PP_REPEAT_FROM_TO(1,6,GENERATE_TUPLE_SERIALIZE,~);
#   undef GENERATE_TUPLE_SERIALIZE

    template< typename Archive, class T1, class T2 > 
    void serialize(Archive & ar,                                        
                   boost::tuples::cons<T1, T2> & t,   
                   const unsigned int version)                         
    {                                                                 
      ar & boost::tuples::get<0>(t);                                  
      ar & boost::tuples::get<1>(t);                                 
    }


    template< typename Archive, typename T > 
    void serialize(Archive & ar,                                            
                   boost::tuples::cons<T, boost::tuples::null_type > & t,   
                   const unsigned int version)                              
    {                                                                       
      ar & boost::tuples::get<0>(t);                                 
    }

}}

#endif
