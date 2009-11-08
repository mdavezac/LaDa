//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/xpressive/regex_primitives.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_compiler.hpp>

#include <print/manip.h>

#include "mlclusters.h"


namespace LaDa
{
  namespace CE 
  {
    void MLClusters :: init(MLCluster const &_cls, types::t_real _eci)
    {
      eci = _eci;
      resize(1);
      back() = _cls;

      // checks for existence of lattice.
      LADA_DOASSERT(Crystal::Structure::lattice,"No lattice found.\n");
      Crystal::Lattice lattice = *Crystal::Structure::lattice;
      if( lattice.space_group.size() == 0 ) lattice.find_space_group();
      if( lattice.space_group.size() == 0 ) return;

      Crystal::Lattice::t_SpaceGroup :: const_iterator i_op = lattice.space_group.begin();
      Crystal::Lattice::t_SpaceGroup :: const_iterator const i_op_end = lattice.space_group.end();
      for(; i_op != i_op_end; ++i_op)
      {
        MLCluster cluster = _cls; 
        cluster.apply_symmetry(*i_op);
        if( end() == std::find(begin(), end(), cluster) ) push_back(cluster);
      }
    }

    bool MLClusters :: Load( const TiXmlElement &_node )
    {
      clear();
      if( not _node.Attribute("eci") ) eci = 0e0;
      else eci = boost::lexical_cast<types::t_real>(_node.Attribute("eci"));

      const TiXmlElement *child = _node.FirstChildElement( "Cluster" );
      MLCluster cluster;
      for ( ; child; child=child->NextSiblingElement( "Cluster" ) )
      {
        if( not cluster.Load(*child) ) return false;
        push_back(cluster);
      }

      return true;
    }

    bool MLClusters::operator==(MLClusters const &_b) const
    {
      if( size() != _b.size() ) return false;
      MLClusters::const_iterator i_cls = _b.begin();
      MLClusters::const_iterator const i_cls_end = _b.end();
      for(; i_cls != i_cls_end; ++i_cls)
        if( end() == std::find(begin(), end(), *i_cls) ) return false;
      return true;
    }
    
    types::t_real MLClusters::operator()( Crystal::Structure const &_str, 
                                          std::vector< std::vector<size_t> > const &_map,
                                          Crystal::t_SmithTransform const &_transform ) const
    {
      if( empty() ) return 1;

      types::t_real result(0);
      const_iterator i_cluster = begin();
      const_iterator const i_cluster_end = end();
      for(; i_cluster != i_cluster_end; ++i_cluster) // loop over clusters classes.
        result += (*i_cluster)(_str, _map, _transform);
      return result;
    }

    std::ostream& operator<<( std::ostream &_stream,  MLClusters const &_class )
    {
      if( _class.size() == 0 )
        return _stream << "Cluster class J0: " << _class.eci << "\n";

      _stream << "Cluster Class: " << _class.eci << "\n";
      std::ostringstream stream;
      MLClusters::const_iterator i_cls = _class.begin();
      MLClusters::const_iterator const i_cls_end = _class.end();
      for(; i_cls != i_cls_end; ++i_cls) stream << *i_cls;
      return _stream << Print::indent( stream.str(), "  ");
    }
    
    bool bypass_comment( std::istream & _sstr, std::string &_line )
    {
      namespace bx = boost::xpressive;
      bx::sregex const re = bx::bos >> (*bx::_s) >> bx::as_xpr('#');
      do { std::getline( _sstr, _line ); }
      while(     bx::regex_search( _line, re ) 
             and ( not _sstr.eof() ) );
      return not _sstr.eof();
    };

    bool read_cluster( const Crystal::Lattice &_lat, 
                       std::istream & _sstr,
                       MLCluster &_out )
    {
      _out.clear(); // clears spins.

      std::string line;
      // bypass comments until a name is found.
      do
      { // check for comment.
        if( not bypass_comment(_sstr, line) ) return false;
      } while( not (    line.find( 'B' ) != std::string::npos 
                     or line.find( 'J' ) != std::string::npos  ) );

      // now read cluster.
      if( not bypass_comment(_sstr, line) ) return false;
      types::t_unsigned order;
      types::t_unsigned multiplicity;
      { // should contain the order and the multiplicity.
        std::istringstream sstr( line ); 
        sstr >> order >> multiplicity;
        LADA_DOASSERT( not sstr.bad(), "Error while reading figure.\n" );
      }

      // creates cluster.
      _out.origin.site = 0;
      _out.origin.pos = _lat.sites[0].pos;
      // bypasse first vector (eg origin).
      LADA_ASSERT( bypass_comment(_sstr, line), "Unexpected end of file.\n");
      for(; order; --order)
      {
        LADA_ASSERT( bypass_comment(_sstr, line), "Unexpected end of file.\n");
        std::istringstream sstr(line);
        MLCluster::Spin spin;
        spin.site = 0;
        sstr >> spin.pos(0) >> spin.pos(1) >> spin.pos(2);
        _out.push_back( spin );
      }
      return true;
    }

    boost::shared_ptr<t_MLClusterClasses> read_clusters( const Crystal::Lattice &_lat, 
                                                         const boost::filesystem::path &_path, 
                                                         const std::string & _genes ) 
    {
      namespace fs = boost::filesystem;  
      // Check path.
      LADA_DOASSERT( fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      LADA_DOASSERT( fs::is_regular( _path ) or fs::is_symlink( _path ),
                     _path << " is neither a regulare file nor a system link.\n" )

      // create result.
      boost::shared_ptr<t_MLClusterClasses> result( new t_MLClusterClasses );

      std::ifstream file( _path.string().c_str(), std::ifstream::in );
      std::string line;

      // create typical cluster.
      MLCluster cluster;

      // read new clusters until done.
      size_t i(0);
      while( read_cluster( _lat, file, cluster ) )
      {
        // This is a hack to avoid inputing Zhe's tetragonal pair terms.
        if( cluster.order() == 2 ) continue;

        ++i;
        if( not _genes.empty() )
        {
          LADA_DOASSERT( i <= _genes.size(), "Gene size and jtypes are inconsistent.\n" )
          if( _genes[i-1] == '0' ) continue; 
        }

        result->push_back(MLClusters());
        result->back().init(cluster);
      }

      LADA_DOASSERT( _genes.empty() or i == _genes.size(), "Gene size and jtypes are inconsistent.\n" )
    }


  } // namespace CE
}
