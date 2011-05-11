//
//  Version: $Id: control_result_section.impl.h 1227 2009-07-14 02:17:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_CONTROL_RESULT_SECTION_H_
# define _LADA_LOADNSAVE_INITIALIZER_CONTROL_RESULT_SECTION_H_

  struct FlowControl
  {
    //! Whether the section was found in xml tree
    bool is_in_xml_tree;
    //! Whether the section is required.
    bool is_required;
    //! Whether found all required options.
    bool found_required_options;
    //! Whether found all required subsections.
    bool found_required_subsections;
    //! Whether found all required subsections.
    bool parsed_options;
    //! Whether found all required subsections.
    bool parsed_subsections;
    //! Whether found all required subsections.
    bool found;
    //! Does not print errors.
    bool quiet;
    //! Does not check for unknown sections
    bool check_unknown_sections;
    //! Does not check for unknown options
    bool check_unknown_options;

    //! Default constructor.
    FlowControl() : is_in_xml_tree(false),
                    is_required(false),
                    found_required_options(true),
                    found_required_subsections(true),
                    parsed_subsections(true),
                    parsed_options(true),
                    found(false),
                    quiet(false),
                    check_unknown_sections(true),
                    check_unknown_options(true) {}

    //! Copy constructor.
    FlowControl   ( FlowControl const &_c )
                : is_in_xml_tree(_c.is_in_xml_tree),
                  is_required(_c.is_required),
                  found_required_options(_c.found_required_options),
                  found_required_subsections(_c.found_required_subsections),
                  parsed_subsections(_c.parsed_subsections),
                  parsed_options(_c.parsed_options),
                  found(_c.found),
                  quiet(_c.quiet),
                  check_unknown_sections(_c.check_unknown_sections),
                  check_unknown_options(_c.check_unknown_options) {}

    //! Prints error if any, returns true or false if parsed.
    template< class T_EXPR > bool result( const T_EXPR& _expr ) const;
    //! Whether to parse section.
    template< class T_EXPR > bool do_parse( const T_EXPR& _expr, tree :: Section const& _sec );
    //! Whether parsing went ok.
    void did_parse(std::string const& _name);
  };
#elif !defined(_LADA_LOADNSAVE_INITIALIZER_CONTROL_RESULT_SECTION_BODY_H_)
# define _LADA_LOADNSAVE_INITIALIZER_CONTROL_RESULT_SECTION_BODY_H_

  template< class T_EXPR >
    bool Section_::FlowControl::result( const T_EXPR& _expr ) const
    {
      std::string const name( transform::GetName()(_expr) );

      typedef typename boost::result_of< transform::GetContent(T_EXPR) >::type t_Content;
      t_Content const content( transform::GetContent()(_expr) );
      TaggedOptions req_options( content, tags::option::required );
      TaggedSections req_sections( content, tags::section::required );
      TaggedOptions const idops( content, tags::option::id );

      if( not found )
      {
        if( not is_required ) return true;
        if( quiet ) return false;
        if( not is_in_xml_tree )
             std::cerr << "Required section " << name << " not found.\n";
        else std::cerr << "Could not find a " << name
                       << " section with all id options,"
                          " required options, and required sub-sections:\n"
                       << "  Id options: " << idops 
                       << "  Required options: " << req_options 
                       << "  Required subsections: " << req_sections;
        return false;
      }

      if( not quiet )
      {
        if( not found_required_options ) 
          std::cerr << "Could not find all required options in " 
                    << name << ": " << req_options;
        if( not found_required_subsections ) 
          std::cerr << "Could not find all required subsections in " 
                    << name << ": " << req_sections;
        if( not parsed_options )
          std::cerr << "Could not parse options in " + name + ".\n";
        if( not parsed_subsections )
          std::cerr << "Could not parse subsections in " + name + ".\n";
      }
      return     found_required_options 
             and found_required_subsections
             and parsed_options
             and parsed_subsections;
    };

  template< class T_EXPR >
    bool Section_::FlowControl :: do_parse( const T_EXPR& _expr,
                                            tree :: Section const& _sec ) 
    {
      std::string const name( transform::GetName()(_expr) );
      typedef typename boost::result_of< transform::GetContent(T_EXPR) >::type t_Content;
      t_Content const content( transform::GetContent()(_expr) );
      TaggedOptions req_options( content, tags::option::required );
      TaggedSections req_sections( content, tags::section::required );

      if( not req_options( _sec ) )  found_required_options = false;
      if( not req_sections( _sec ) ) found_required_subsections = false;

      return found_required_subsections and found_required_options;
    }


#endif
