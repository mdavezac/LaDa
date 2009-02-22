//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include "../functional.h"
#include "../layered.h"
#include "../va.h"

#include "vff.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      template<class T_TYPE>
        void Vff_from_XML(T_TYPE &_type, const std::string &_filename )
        {
          std::cout << "lattice: " << *_type.structure.lattice << "\n";
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
        
          __DOASSERT( not doc.LoadFile(), 
                         doc.ErrorDesc() << "\n"  
                      << "Could not load input file " << _filename  
                      << ".\nAborting.\n" ) 
          __DOASSERT( not docHandle.FirstChild("Job").Element(),
                      "Could not find <Job> tag in " << _filename << ".\n" )
         
          std::cout << "type: " << _type.structure << "\n";
          __DOASSERT( not _type.Load( *docHandle.FirstChild("Job").Element() ),
                         "Could not load Vff functional from " + _filename + ".\n" )
        }
    }


    //! Assumes ownership of the Crystal::Structure object needed by vff.
    template< class T_VFF > class Vff 
    {
      public:
        //! Type of the functional.
        typedef LaDa::Vff::VABase< T_VFF > t_Functional;
        //! Constructor.
        Vff() { functional_.reset( new t_Functional( structure ) ); }
        //! Copy Constructor.
        Vff( const Vff& _c ) : structure( _c.structure )
          { functional_.reset( new t_Functional( structure ) ); }

        //! Redo initialization of vff from scratch.
        void init() { functional_->init( true ); }

        //! fast evaluation with no reinitialization.
        types::t_real operator()() const { return functional_->evaluate() / 16.0217733; }
        //! Returns the stress.
        atat::rMatrix3d get_stress() const { return functional_->get_stress(); }

        //! Loads from an XML input file.
        bool Load( const TiXmlElement &_node ) { return functional_->Load( _node ); }

        //! Prints escan input.
        void print_escan_input( const std::string &_path ) const
          { functional_->print_escan_input( _path ); }
          
        //! The owned structure.
        Crystal::Structure structure;


      protected:
        //! The functional.
        boost::shared_ptr< LaDa::Vff::VABase< T_VFF > > functional_;
    };

    class LVff : public LaDa::Vff::Layered
    {
      public:
        //! Constructor.
        LVff( LaDa::Crystal::Structure& _str ) : LaDa::Vff::Layered( _str ) {}
        //! Exposes direction setters.
        atat::rVector3d get_direction() const { return direction; } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
        {
          is_fixed_by_input = true;
          direction = _direction;
          create_template_strain();
        } 
    };

    class LayeredVff : public Vff< LVff >
    {
      public:
        //! Exposes direction setters.
        atat::rVector3d get_direction() const
          { return functional_->Vff().get_direction(); } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
          { functional_->Vff().set_direction( _direction ); } 
    };


#   ifdef EXPOSEVFF 
#     error Macro EXPOSEVFF already exists.
#   endif 
#   define EXPOSEVFF( a, b, c ) \
      bp::class_< b >( a, c ) \
        .def( bp::init< b >() ) \
        .def_readwrite( "structure",    &b::structure, \
                        "The structure to minimize. Input and Output to the functional." ) \
        .def( "fromXML",  &XML::Vff_from_XML<b>, bp::arg("file"),\
              "Loads the vff parameters from an XML file." ) \
        .def( "evaluate",  &b::operator(), \
              "Minimizes the current structure and returns the energy in eV." ) \
        .def( "init",  &b::init, "Initializes the functional for the current structure." ) \
        .def( "print_escan_input",  &b::print_escan_input, bp::arg("file"), \
              "Outputs the current structure in a format suitable for pescan." ) \
        .add_property( "stress",  &b::get_stress, \
                       "Returns the stress. Meaningfull only following a call to Vff.evaluate()." ) 

    void expose_vff()
    {
      typedef Vff< LaDa::Vff::Functional > t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF
      ( 
        "Vff", 
        t_Vff, 
        "A Valence Force Field Functional.\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :)."
      );
    }

    void expose_layeredvff()
    {
      typedef LayeredVff t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF
      ( 
        "LayeredVff", 
        t_Vff,
        "A Valence Force Field Functional with epitaxial constraints.\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :)."
      ).add_property( "direction",  &t_Vff::get_direction, &t_Vff::set_direction,
                      "Defines the direction of growth." );
    }

#   undef EXPOSEVFF

  }
} // namespace LaDa
