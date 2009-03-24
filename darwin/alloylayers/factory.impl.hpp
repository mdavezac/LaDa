//
//  Version: $Id$
//

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
      //! Holds some functions to connect quantities, functionals, printing policies.
      namespace Factory
      {
        namespace details
        {
          //! Dumps object to a string.
          std::ostream& printg( std::ostream& _stream,
                                const LaDa::GA::BitString::Object<>& _o );
          //! Actually function which does the connecting.
          template< class T_EVALUATOR, class T_PROPERTY >
            class ConnectProperty
            {
              public:
                //! Constructor.
                ConnectProperty   ( T_EVALUATOR& _evaluator,
                                    T_PROPERTY _property, 
                                    const std::string &_name,
                                    const std::string &_decl )
                                : evaluator_( _evaluator ),
                                  property_( _property ), 
                                  name_( _name ), decl_( _decl ) {}
                //! Copy Constructor.
                ConnectProperty   ( const ConnectProperty& _c )
                                : evaluator_( _c.evaluator_ ), 
                                  property_( _c.property_ ), 
                                  name_( _c.name_ ), decl_( _c.decl_ ) {}
  
                //! Functor.
                void operator()();

              private:
                //! Evaluator for connect to.
                T_EVALUATOR &evaluator_;
                //!  Property to connect.
                T_PROPERTY property_;
                //! Name of the property.
                const std::string name_;
                //! Output declaration.
                const std::string decl_;
            };
        }
        template< class T_EVALUATOR, class T_PROPERTY >
          details::ConnectProperty<T_EVALUATOR, T_PROPERTY> 
            connect_property( T_EVALUATOR &_evaluator, T_PROPERTY _property, 
                              const std::string &_name, const std::string &_decl )
            {
              return details::ConnectProperty<T_EVALUATOR, T_PROPERTY>
                                             ( _evaluator, _property, _name, _decl );
            }

        template< class T_EVALUATOR > void genotype( const T_EVALUATOR& )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          t_Object :: connect_print( bl::bind( &details::printg, bl::_1, bl::_2 ) ); 
        }

        template< class T_EVALUATOR > void strain_energy( T_EVALUATOR & _evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            ( 
              &t_Quantity::push_back,
              bl::_2,
              bl::bind( &t_Object::energy, bl::_1 ) / bl::constant(16.0217733) 
            )
          );
          t_Object :: connect_print
          (
            bl::_1 << bl::constant(" Strain Energy: ")
                   << bl::bind( &t_Object::energy, bl::_2 ) / bl::constant( 16.0217733 )
                   << bl::constant( " eV/f.u. " )
          );
          declare( "Strain energy(VFF)" );
        }

        template< class T_EVALUATOR > void epitaxial_strain( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
              bl::_2,
              bl::bind
              ( 
                &Vff::inplane_stress,
                bl::bind<atat::rMatrix3d>( &t_Object :: stress, bl::_1 ),
                bl::bind<atat::rVector3d>( &T_EVALUATOR :: get_direction,
                                           bl::constant( boost::cref( _evaluator ) ) )
              ) / bl::constant( 16.0217733 )
            )
          );
          t_Object :: connect_print
          (
            bl::_1 << bl::constant(" Epi Strain: ")
                   << bl::constant( std::fixed ) << bl::constant( std::setw(8) )
                   << bl::constant( std::setprecision(3) )
                   << bl::bind
                      ( 
                        &Vff::inplane_stress,
                        bl::bind<atat::rMatrix3d>( &t_Object :: stress, bl::_2 ),
                        bl::bind<atat::rVector3d>( &T_EVALUATOR :: get_direction,
                                                   bl::constant( boost::cref( _evaluator ) ) )
                      ) / bl::constant( 16.0217733 ) << bl::constant( " eV/f.u. " )
          );
          declare( "Epitaxial strain (VFF)" );
        }

        template< class T_EVALUATOR > void bandgap( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
             bl::_2,
              bl::bind( &t_Object::vbm, bl::_1 ) - bl::bind( &t_Object::cbm, bl::_1 )
            )
          );
          t_Object :: connect_print
          ( 
            bl::_1 << bl::constant( "Bandgap: ") 
                   << bl::bind( &t_Object::vbm, bl::_2 ) - bl::bind( &t_Object::cbm, bl::_2 )
                   << bl::constant( " " )
          );
          declare( "Band-gaps(Escan)" );
        }

        template< class T_EVALUATOR > void dipole( T_EVALUATOR &_evaluator )
        {
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;

          _evaluator.do_dipole = true;
          // connect_property returns a functor. Next line calls the functor.
          connect_property( _evaluator, &t_Object::osc_strength,
                            "transition dipole=", 
                            "Dipole oscillator strength between "
                            "VBM and CBM. Arbitrary units." )();
        }

        template< class T_EVALUATOR > void emass( T_EVALUATOR &_evaluator )
        {
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;

          _evaluator.do_emass = true;
          // connect_property returns a functor. Next line calls the functor.
          connect_property( _evaluator, &t_Object::emass,
                            "emass=", 
                            "Electron effective mass (a.u.)" )();
        }
        template< class T_EVALUATOR > void hmass( T_EVALUATOR &_evaluator )
        {
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;

          _evaluator.do_hmass = true;
          connect_property( _evaluator, &t_Object::hmass,
                            "hmass=", 
                            "Hole effective mass (a.u.)" )();
        }

        void declare( const std::string& _string )
        {
          LaDa::Print::out << _string << " will be optimized.\n";
          LaDa::Print::xmg << LaDa::Print::Xmg::comment
                           <<  _string << " will be optimized.\n";
        }

        template< class T_FACTORY >
          void read_physical_properties( T_FACTORY& _factory,
                                         boost::filesystem::path &_path )
          {
            TiXmlDocument doc; 
            std::string filetxt;
            opt::read_xmlfile( _path, filetxt );
            doc.Parse( filetxt.c_str() );
            TiXmlHandle docHandle( &doc ); 

            const TiXmlElement* child = docHandle.FirstChild( "Job" )
                                                 .FirstChild( "GA" )
                                                 .FirstChild( "Objectives" )
                                                 .FirstChild( "Objective" )
                                                 .Element();
            __DOASSERT( not child, "Could not find Objective tag.\n" )
            for(; child; child = child->NextSiblingElement("Objective") )
            {
              __DOASSERT( not child->Attribute("value"),
                          "Found objective tag without a value attribute.\n" )
              const std::string value = child->Attribute( "value" );
              _factory( value );
            }
          }

        namespace details
        {
          template< class T_EVALUATOR, class T_PROPERTY >
            void ConnectProperty<T_EVALUATOR, T_PROPERTY > :: operator()()
            {
              namespace bl = boost::lambda;
              typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
              typedef typename t_IndivTraits :: t_Object t_Object;
              typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
              evaluator_.connect
              (
                bl::bind
                (
                  &t_Quantity::push_back,
                  bl::_2,
                  bl::bind( property_, bl::_1 )
                )
              );
              t_Object :: connect_print
              ( 
                bl::_1 << bl::constant(name_) << bl::bind( property_, bl::_2 )
                       << bl::constant( " " )
              );
              declare( decl_ );
            }

          inline std::ostream& printg( std::ostream& _stream,
                                       const LaDa::GA::BitString::Object<>& _o )
          {
            foreach( types::t_real var, _o.Container() )
                _stream << std::fixed << std::setw(7) << std::setprecision(3) << var;
            return _stream;
          }
        }

      } // namespace factory
    } // namespace Layered

  }
} // namespace LaDa

