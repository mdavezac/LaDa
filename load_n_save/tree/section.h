#ifndef LADA_LNS_PARSED_TREE_SECTION_H
#define LADA_LNS_PARSED_TREE_SECTION_H

#include "LaDaConfig.h"

#include <map>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/lexical_cast.hpp>

#include "../string_type.h"
#include "base.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace tree
    {
      //! An option holds a name and a value.
      struct Option
      {
        //! Name of the option.
        t_String name;
        //! Value of the option.
        t_String value;
        //! Wether this option has been parsed yet,
        mutable bool fresh;

        //! Constructor.
        Option() : name(""), value(""), fresh(true) {}
        //! Copy Constructor.
        Option   ( Option const& _c )
               : name(_c.name), value(_c.value), fresh(_c.fresh) {}
      };

      
      //! \brief Xml node.
      //! \details Holds xml tag name, options and subnodes. Allows parsing from
      //!          string and direct manipulation of the xml tree.
      class Section : public virtual details::Base<Section>, public virtual details::Options<Option>
      {
          //! Type of option tree.
          typedef details::Options<Option> t_Options;
          //! Type of option tree.
          typedef details::Base<Section> t_Subsections;
        public:
          //! Iterators
          struct iterator
          {
            //! Type of the subsections iterator.
            typedef t_Subsections :: iterator :: subsection subsection;
            //! Type of the option iterator.
            typedef t_Options :: iterator option;
          };
          //! Iterators
          struct const_iterator
          {
            //! Type of the subsections iterator.
            typedef t_Subsections :: const_iterator :: subsection subsection;
            //! Type of the subsections iterator.
            typedef t_Options :: const_iterator option;
          };
    
          //! Tag Name.
          t_String name;
          //! Whether this section has not been read yet.
          mutable bool fresh;

          //! Constructor.
          Section( const t_String _name ) : name( _name ), fresh(true) {}
          //! Copy Constructor.
          Section   ( Section const& _c )
                  : t_Subsections(_c), t_Options(_c), name(_c.name), fresh(_c.fresh) {}
          //! Copy Constructor.
          Section   ( t_Subsections const& _c )
                  : t_Subsections(_c), t_Options(), name(""), fresh(true) {}


          //! Insert an option.
          void push_back( t_String const& _name, t_String const &_value = "" )
          { 
            Option op; op.name = _name; op.value = _value;
            t_Options::push_back( op);
          }
          //! Insert an option.
          template< class T_TYPE >
            typename boost::disable_if< typename boost::is_same< T_TYPE, t_String > :: type > :: type
              push_back( t_String const& _name, T_TYPE const &_value )
              {
                t_String const value = boost::lexical_cast<t_String>( _value );
                push_back( _name, value );
              }
          
          //! Copies a section.
          void copy_to( Section& _c ) const
            { t_Subsections::copy_to(_c); t_Options::copy_to(_c); _c.name = name; _c.fresh = fresh; }
      };

      typedef details::Base<Section> Base;

    } // namespace parser

  } // namespace load_n_save
} // namespace LaDa

#endif
