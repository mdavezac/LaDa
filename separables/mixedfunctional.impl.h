#include <boost/bind.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <iostream>

#include <ce/find_pis.h>
#include <ce/create_pairs.h>

namespace LaDa
{
  namespace CE
  {
#   if defined( HEAD ) || defined(INHEAD) || defined(INHEAD2)
#     error "Macros with same names."
#   endif
#   define HEAD \
      MixedSeparables<T_TRAITS, T_HARMONICS> 
#    define INHEAD( code ) \
       template< class T_TRAITS, class T_HARMONICS >  code HEAD
#    define INHEAD2( code1, code2 ) \
       template< class T_TRAITS, class T_HARMONICS >  code1, code2 HEAD
    INHEAD( template< class TT_TRAITS > void )
      :: operator=( const Separables<TT_TRAITS> &_c )
      {
        t_SepBase::norms = _c.norms;
        t_SepBase::coefficients() = _c.coefficients(); 
      }
    INHEAD( template< class T_COLTRAITS > void )
      :: operator=( const MixedApproach<T_COLTRAITS> &_col )
      {
        namespace bblas = boost::numeric::ublas;
        this->coefficients_interface()
          = _col.separables().coefficients_interface();
        this->operator=( _col.separables() );
        typedef typename MixedApproach<T_COLTRAITS> :: t_Matrix t_MixedAMat;
        const bblas::matrix_column< const t_MixedAMat > column0( _col.coefficients(), 0 );
        ecis().resize( _col.dof() - _col.separables().dof() );
        std::copy( column0.begin() + _col.separables().dof(), column0.end(), ecis().begin() );
        clusterclasses().clear();
        std::copy( _col.clusters().begin(), _col.clusters().end(), 
                   std::back_inserter( clusterclasses() ) );
      }
 
 
    INHEAD( typename HEAD :: t_SepBase :: t_Matrix :: value_type )
      ::  operator()( const Crystal::Structure &_structure )
      {
        namespace bl = boost::lambda;
        namespace bblas = boost::numeric::ublas;
        LADA_DEBUG_TRY_BEGIN
        typename t_SepBase::t_Matrix::value_type result(0);
        typedef t_PosBase :: t_Configurations t_Configurations;
        if( t_PosBase::dof() )
        {
          t_Configurations configs;
          t_PosBase::operator()( _structure, configs );
          t_Configurations :: const_iterator i_config = configs.begin();
          t_Configurations :: const_iterator i_config_end = configs.end();
          for(; i_config != i_config_end; ++i_config )
            result +=   t_SepBase::operator()( i_config->first )
                      * i_config->second;
        }
        if( ecis().size() and ( not clusterclasses().empty() ) )
        {
          typename t_SepBase::t_Vector pis( clusterclasses().size() );
          std::fill( pis.begin(), pis.end(),
                     typename t_SepBase::t_Vector::value_type(0) );
          find_pis( clusterclasses(), _structure, pis );
          result += bblas::inner_prod( pis, ecis() );
        }
        if( t_CSBase::harmonics.size() )
        {
          t_CSBase::operator<<( _structure );
          t_CSBase::resize( _structure.atoms.size() );
          std::transform( _structure.atoms.begin(), _structure.atoms.end(), t_CSBase::begin(),
                          boost::bind( &Crystal::Structure::t_Atom::type, _1 ) );
          result += t_CSBase::evaluate();
        }
        return result;
        LADA_DEBUG_TRY_END(, "Error while evaluating mixed-approach functional.\n" )
      }
 
    INHEAD( void ) :: init( const std::string &_csxml,
                            const std::string &_desc, 
                            const types::t_unsigned _nbpairs,
                            const std::string &_jtypes,
                            bool _rmpairs, bool _addJ0, bool _addJ1 )
    {
      LADA_DEBUG_TRY_BEGIN
      const std::string csxml  = boost::algorithm::trim_copy( _csxml );
      const std::string desc   = boost::algorithm::trim_copy( _desc );
      const std::string jtypes = boost::algorithm::trim_copy( _jtypes );
      if( not csxml.empty() )
      {
        TiXmlDocument doc( _csxml );
        LADA_DO_NASSERT( not doc.LoadFile(), 
                      "error while opening input file "
                   << _csxml << "\n" << doc.ErrorDesc()  )
        TiXmlHandle handle( &doc );
        const TiXmlElement *child = handle.FirstChild( "Job" ).FirstChild( "CS" ).Element();
        LADA_DO_NASSERT( not child, "Could not find CS in input.\n" )
        LADA_DO_NASSERT( not t_CSBase::Load_Harmonics( *child ),
                    "Error while loading harmonics from input.\n" )
      }
      if( not desc.empty() ) t_PosBase :: create_positions( desc ); 
      std::vector< t_ClusterClass > cl;
      if( not _nbpairs == 0 )
        create_pairs( *Crystal::Structure::lattice, _nbpairs, cl );
      if( not jtypes.empty() )
        CE::read_clusters( *Crystal::Structure::lattice, jtypes, cl ); 
      clusterclasses().clear();
      if( not cl.size() == 0 )
        std::copy( cl.begin(), cl.end(), std::back_inserter( clusterclasses() ) );
      CE::Cluster cluster;
      if( _addJ0 )
      {
        bool isfound = false;
        foreach( t_ClusterClass &_class, clusterclasses() )
          if( _class.front().size() == 0 ) { isfound = true; break; }
        if( isfound ) clusterclasses().push_back( std::vector<CE::Cluster>(1, cluster) ); 
      }
      cluster.Vectors().resize(1, math::rVector3d(0,0,0) );
      if( _addJ1 )
      {
        bool isfound = false;
        foreach( t_ClusterClass &_class, clusterclasses() )
          if( _class.front().size() == 1 ) { isfound = true; break; }
        if( isfound ) clusterclasses().push_back( std::vector<CE::Cluster>(1, cluster) ); 
      }
      ecis().resize( clusterclasses().size() );
      LADA_DEBUG_TRY_END(, "Error in MixedFunctional::init().\n" )
    }
 
    template<class T_TRAITS, class T_HARMONICS> 
      std::ostream& operator<<( std::ostream& _stream, const HEAD &_sep )
      {
        typedef typename HEAD :: t_SepBase  t_SepBase;
        typedef typename HEAD :: t_CSBase   t_CSBase;
        typedef typename HEAD :: t_PosBase  t_PosBase;
        std::cout << "Mixed Functional:\n";
        if( _sep.t_PosBase::dof() ) _stream << *( (t_SepBase*) &_sep ) << "\n";
        if( _sep.ecis().size() and ( not _sep.clusterclasses().empty() ) )
          foreach( const typename HEAD :: t_ClusterClass &_class, _sep.clusterclasses() )
            _stream << _class.front()  << " D=" << _class.size() << "\n";
        if( _sep.t_CSBase::harmonics.size() ) _stream << " CS. no description. \n";
        return _stream;
      }
#    undef HEAD
#    undef INHEAD
#    undef INHEAD2
  } // end of CE namespace.
} // namespace LaDa
