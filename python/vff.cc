//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <vff/functional.h>
#include <vff/layered.h>
#include <vff/va.h>

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
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
        
          __DOASSERT( not doc.LoadFile(), 
                         doc.ErrorDesc() << "\n"  
                      << "Could not load input file " << _filename  
                      << ".\nAborting.\n" ) 
          __DOASSERT( not docHandle.FirstChild("Job").Element(),
                      "Could not find <Job> tag in " << _filename << ".\n" )
         
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
        types::t_real operator()() const { return functional_->evaluate(); }
        //! Returns the stress.
        atat::rMatrix3d get_stress() const { return functional_->get_stress(); }

        //! Loads from an XML input file.
        bool Load( const TiXmlElement &_node ) { return functional_->Load( _node ); }
       
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
        atat::rVector3d get_direction() const { return functional_->Vff().get_direction(); } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
          { functional_->Vff().set_direction( _direction ); } 
    };


#   ifdef EXPOSEVFF 
#     error Macro EXPOSEVFF already exists.
#   endif 
#   define EXPOSEVFF( a, b ) \
      bp::class_< b >( a ) \
        .def( bp::init< b >() ) \
        .def_readwrite( "structure",    &b::structure ) \
        .def( "fromXML",  &XML::Vff_from_XML<b> ) \
        .def( "evaluate",  &b::operator() ) \
        .def( "init",  &b::init ) \
        .add_property( "stress",  &b::get_stress ) 

    void expose_vff()
    {
      typedef Vff< LaDa::Vff::Functional > t_Vff;
      namespace bp = boost::python;
     //bp::class_< t_Vff >( "Vff" ) 
     //  .def( bp::init< t_Vff >() )
     //  .def_readwrite( "structure", &t_Vff::structure )
     //  .def( "fromXML",  &XML::Vff_from_XML<t_Vff> ) 
     //  .def( "evaluate",  &t_Vff::operator() ) 
     //  .def( "init",  &t_Vff::init )
     //  .add_property( "stress",  &t_Vff::get_stress );
      EXPOSEVFF( "Vff", t_Vff );
    }

    void expose_layeredvff()
    {
      typedef LayeredVff t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF( "LayeredVff", t_Vff )
        .add_property( "direction",  &t_Vff::get_direction, &t_Vff::set_direction );
    }

#   undef EXPOSEVFF

  }
} // namespace LaDa
