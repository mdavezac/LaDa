//
//  Version: $Id: option_.cc 1226 2009-07-13 06:28:01Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "option_.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      Option_::FlowControl :: FlowControl() : match_only(false),
                                              match_tag( tags::option::default_ ),
                                              is_many(false),
                                              is_id(false),
                                              is_required(false),
                                              found(false), 
                                              found_many(false),
                                              is_in_xml_tree(false),
                                              no_default(false) {};
      Option_::FlowControl :: FlowControl   ( FlowControl const& _c )
                                          : match_only(_c.match_only),
                                            match_tag( _c.match_tag ),
                                            is_many(_c.is_many),
                                            is_id(_c.is_id),
                                            is_required(_c.is_required),
                                            found(_c.found),
                                            found_many(_c.found_many),
                                            is_in_xml_tree( _c.is_in_xml_tree ),
                                            no_default(_c.no_default) {};
      //! True if matches correct tag
      bool Option_::FlowControl :: good_tag( size_t _tag ) const
      {
        if( match_tag == tags::option::default_ ) return true;
        return _tag & match_tag;
      }

//     inline void Option_::FlowControl :: found_another() const
//     {
//       if( not found ) found = true;
//       else found_many = true;
//     }

      void Option_::FlowControl :: init( size_t _tag )
      {
        found = false;
        found_many = false;
        is_many = _tag & tags::option::many;
        is_id = _tag & tags::option::id;
        is_required = _tag & tags::option::required;
      }

      bool Option_::FlowControl ::  operator()( std::string const& _name ) const
      {
        if( match_only ) return found;

        if( found )
        {
          if( not found_many ) return true;
          if( is_many ) return true;
          std::cerr << "Found more than one option " + _name +
                       "When only one is expected.\n";
          return false;
        }

        if( is_id ) return false;


        if(is_in_xml_tree)
          std::cerr << "Error while parsing option " + _name + ".\n";
        if( is_required)
          std::cerr << "Could not find/parse required option " + _name + ".\n";

        return not(is_required or is_in_xml_tree);
      };
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


