//
//  Version: $Id$
//

#include <boost/lambda/bind.hpp>
#include <iostream>

namespace CE
{
  template< class T_TRAITS > template< class TT_TRAITS >
    void MixedSeparables<T_TRAITS> :: operator=( cont Separables<TT_TRAITS> &_c )
    {
      t_SepBase::norms = _c.norms;
      t_SepBase::coefficients() = _c.coefficients(); 
    }

  template< class T_TRAITS > 
    typename MixedSeparables<T_TRAITS>::t_SepBase::t_Matrix::value_type
      MixedSeparables<T_TRAITS> ::  operator()( const Crystal::Structure &_structure )
      {
        namespace bl = boost::lambda;
        __DEBUGTRYBEGIN
        typename t_SepBase::t_Matrix::value_type result(0);
        typedef t_PosToConf :: t_Configurations t_Configurations;
        if( t_PosBase::positions.size() )
        {
          t_Configurations configs;
          t_PosBase::operator()( _structure, configs );
          t_Configurations :: const_iterator i_config = configs.begin();
          t_Configurations :: const_iterator i_config_end = configs.end();
          for(; i_config != i_config_end; ++i_config )
            result +=   t_SepBase::operator()( i_config->first )
                      * i_config->second;
        }
        if( ecis.size() and ( not clusterclasses().empty() ) )
        {
          typename t_SepBase::t_Vector pis( clusterclasses().size() );
          std::fill( pis.begin(), pis.end(),
                     typename t_SepBase::t_Vector::value_type(0) );
          find_pis( clusterclasses(), _structure, pis );
          result += bblas::inner_prod( pis, ecis() );
        }
        if( t_CSBase::harmonics->size() )
        {
          t_CSBase::operator<<( _structure );
          result += t_CSBase::operator()( _structure.get_concentration() );
        }
        return result;
        __DEBUGTRYEND(, "Error while evaluating mixed-approach functional.\n" )
      }

  template<class T_SEPTRAITS, class T_COLTRAITS > 
  void operator=( MixedSeparables<T_SEPTRAITS> &_s, 
                  const MixedApproach<T_COLTRAITS> &_col )
  {
    __ASSERT( _col.dof() != _s.dof(), "Inconsistent sizes.\n" )
    __ASSERT( _col.separables().dof() != t_SepBase::dof(), "Inconsistent sizes.\n" )
    _s = _col.separables();
    typedef MixedApproach<T_COLTRAITS> :: t_Matrix t_MixedAMat;
    const bblas::matrix_column< const t_Matrix > column0( _sep.coefficients(), 0 );
    std::copy( column0.begin() + t_SepBase::dof(), column0.end(), ecis.begin() );
  }

  void init( const std::string &_csxml,
             const std::string &_desc, 
             const types::t_unsigned _nbpairs = 0,
             const std::string &_jtypes = "",
             bool _rmpairs = false,
             bool _addJO = false,
             bool _addJ1 = false )
  {
    const std::string csxml = Print::StripEdges( csxml );
    const std::string desc = Print::StripEdges( _desc );
    const std::string jtype = Print::StripEdges( _jtypes );
    if( not csxml.empty() )
    {
      TiXmlDocument doc( filename );
      __DOASSERT( not doc.LoadFile(), 
                    "error while opening input file "
                 << filename << "\n" << doc.ErrorDesc()  )
      TiXmlHandle handle( &doc );
      child = parent->FirstChildElement( "CS" );
      __DOASSERT( not child,
                  "Could not find CS in input.\n" )
      __DOASSERT( not harmonics->Load_Harmonics( *child ),
                  "Error while loading harmonics from input.\n" )
    }
    if( not desc.empty() ) t_PosBase :: create_positions( desc ); 
    if( not jtypes.empty() )
    create_pairs( t_PosBase::lattice, _nbpairs, clustersclasses() );
    CE::read_clusters( *lattice, jtypes, clusterclasses() ); 
    CE::Cluster cluster;
    if( _addJ0 )
    {
      t_ClusterClasses :: const_iterator i_class = clusterclasses().begin();
      t_ClusterClasses :: const_iterator i_class_end = clusterclasses().end();
      for(; i_class != i_class_end; ++i_class )
        if( i_class->front().size() == 0 ) break;
      if( i_class == i_class_end )
        clustersclasses().push_back( std::vector<CE::Cluster>(1, cluster) ); 
    }
    cluster.Vectors().resize(1, atat::rVector3d(0,0,0) );
    if( _addJ1 )
    {
      t_ClusterClasses :: const_iterator i_class = clusterclasses().begin();
      t_ClusterClasses :: const_iterator i_class_end = clusterclasses().end();
      for(; i_class != i_class_end; ++i_class )
        if( i_class->front().size() == 1 ) break;
      if( i_class == i_class_end )
        clustersclasses().push_back( std::vector<CE::Cluster>(1, cluster) ); 
    }
  }

  template< class T_TRAITS, class T_HARMONICS >
    void evaluate_pifile( const std::string &_file,
                          MixedSeparables< T_TRAITS, T_HARMOMICS > &_op )
    {
      Crystal :: Structure;
      std::ifstream file( _file.c_str(), std::ifstream::in );
      do
      {
        if( not Crystal :: read_pifile_structure( file, structure ) ) continue;
        std::cout << _op( structure ) << "\n";
        std::for_each
        (
          structure.atoms.begin(), structure.atoms.end(),
          bl::bind( &Crystal::Structure::t_Atom::pos, bl::_1 ) *= bl::constant(-1e0)
        );
        std::cout << _op( structure ) << "\n";
      }
      while( not file.eof() );
    }
} // end of CE namespace.
