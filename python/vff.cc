//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>

#include <vff/functional.h>
#include <vff/layered.h>

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
         
          __DOASSERT( not _type.Load( not docHandle.FirstChold("Job").Element ),
                         "Could not load " << nodename<T_TYPE>()
                      << " from " << _filename << ".\n" )
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

        //! The owned structure.
        Crystal::Structure structure;

      protected:
        //! The functional.
        boost::shared_ptr< LaDa::Vff::VABase< T_VFF > > functional_;
    };

    //! Exposes direction parameter in a relatively safe way.
    class LayeredVff : public LaDa::Vff::Layered
    {
      public:
        //! Returns direction.
        atat::rVector3d get_direction() const { return direction; } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
        {
          is_fixed_by_input = true;
          direction = _direction;
          create_template_strain();
        } 
    };

#   ifdef EXPOSEVFF 
#     error Macro EXPOSEVFF already exists.
#   endif 
#   define EXPOSEVFF( a, b ) \
      bp::class_< b >( #a ) \
        .def( init< b >() ) \
        .def_readwrite( "structure",    &b::structure ) \
        .def( "fromXML",  &XML::Vff_from_XML ) \
        .def( "evaluate",  &b::operator() ) \
        .def( "init",  &b::init() ) \
        .add_property( "stress",  &b::get_stress() ) 

    void expose_vff()
    {
      typedef Vff< Lada::Vff::Functional > t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF( "Vff", t_Vff );
    }

    void expose_layeredvff()
    {
      typedef Vff< LayeredVff > t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF( "LayeredVff", t_Vff )
        .add_property( "stress",  &b::get_direction(), &b::set_direction );
    }

#   undef EXPOSEVFF

  }
} // namespace LaDa
