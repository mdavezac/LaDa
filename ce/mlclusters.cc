#include "LaDaConfig.h"

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/xpressive/regex_primitives.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_compiler.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <opt/tinyxml.h>

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

    bool MLClusters :: load( const TiXmlElement &_node )
    {
      clear();
      if( not _node.Attribute("eci") ) eci = 0e0;
      else eci = boost::lexical_cast<types::t_real>(_node.Attribute("eci"));

      const TiXmlElement *child = _node.FirstChildElement( "Cluster" );
      MLCluster cluster;
      for ( ; child; child=child->NextSiblingElement( "Cluster" ) )
      {
        if( not cluster.load(*child) ) return false;
        push_back(cluster);
      }
      if( size() == 1 ) init(front(), eci);

      return true;
    }

    bool MLClusters::operator==(MLCluster const &_b) const
    {
      if( order() != _b.order() ) return false;
      return end() != std::find(begin(), end(), _b);
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
      LADA_ASSERT( _str.lattice, "Lattice is not set.\n" );
      LADA_ASSERT( _str.atoms.size() % _str.lattice->sites.size() == 0,
                   "Inconsistent structure and lattice.\n" );
      types::t_real const npersite(_str.atoms.size()/_str.lattice->sites.size());
      types::t_real const factor(1e0/types::t_real(npersite));
      for(; i_cluster != i_cluster_end; ++i_cluster) // loop over clusters classes.
        result += (*i_cluster)(_str, _map, _transform);
      return result * factor;
    }

    std::ostream& operator<<( std::ostream &_stream,  MLClusters const &_class )
    {
      if( _class.size() == 0 )
        return _stream << "Cluster class J0: " << _class.eci << "\n";
      if( _class.order() == 1 )
        return _stream << "Cluster class J1: " << _class.eci << "\n";

      _stream << "Cluster Class: " << _class.eci << "\n";
      std::ostringstream stream;
      MLClusters::const_iterator i_cls = _class.begin();
      MLClusters::const_iterator const i_cls_end = _class.end();
      for(; i_cls != i_cls_end; ++i_cls) stream << *i_cls;
      std::string indented = stream.str();
      boost::algorithm::replace_all(indented, "\n", "\n   ");
      indented.erase(indented.size()-3);
      return _stream << "   " << indented;
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

      LADA_DOASSERT( _genes.empty() or i == _genes.size(),
                     "Gene size and jtypes are inconsistent.\n" )
      return result;
    }

    boost::shared_ptr<t_MLClusterClasses> load_old(TiXmlElement const &_node);
    boost::shared_ptr<t_MLClusterClasses> load_new(TiXmlElement const &_node);
    
    boost::shared_ptr<t_MLClusterClasses> load_mlclusters( boost::filesystem::path const &_path, 
                                                           bool is_newformat )
    {
      namespace fs = boost::filesystem;  
      // Check path.
      LADA_DOASSERT( fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      LADA_DOASSERT( fs::is_regular( _path ) or fs::is_symlink( _path ),
                     _path << " is neither a regulare file nor a system link.\n" )

      std::string text;
      opt::read_xmlfile(_path, text);
      TiXmlDocument doc;
      TiXmlHandle handle( &doc );
      doc.Parse( text.c_str() );
      TiXmlElement const *job_node = handle.FirstChild( "Job" ).Element();
      LADA_DOASSERT( job_node, "Could not find Job tag in input.\n" );
      TiXmlElement const *func_node = opt::find_node(*job_node, "Functional", "type", "CE");
      LADA_DOASSERT( func_node, "Could not find Functional CE in input.\n" );

      return is_newformat ? load_new(*func_node): load_old(*func_node);
    }

    boost::shared_ptr<t_MLClusterClasses> load_old(TiXmlElement const &_node)
    {
      
      TiXmlElement const * clusters_node = opt::find_node(_node, "Clusters");
      LADA_DOASSERT( clusters_node, "Could not find Clusters node in input.\n" );
      TiXmlElement const * cluster_node = clusters_node->FirstChildElement("Cluster");
      LADA_DOASSERT( cluster_node, "Could not find Clusters node in input.\n" );
       
      boost::shared_ptr<t_MLClusterClasses> result(new t_MLClusterClasses);
      MLCluster cluster;
      MLClusters clusters;
      cluster.origin.pos = math::rVector3d(0,0,0);
      cluster.origin.site = 0;
      for(; cluster_node; cluster_node = cluster_node->NextSiblingElement("Cluster") )
      {
        LADA_DOASSERT(cluster_node->Attribute("eci"), "eci attribute not found.\n");
        clusters.eci = boost::lexical_cast<types::t_real>(cluster_node->Attribute("eci"));
    
        TiXmlElement const  *vec_node = cluster_node->FirstChildElement("spin");
        if( not vec_node ) // J0
        {
          clusters.clear();
          result->push_back(clusters);
          continue;
        }
        for(cluster.clear(); vec_node; vec_node = vec_node->NextSiblingElement("spin") )
        {
          LADA_ASSERT(vec_node->Attribute("x"), "Missing x attribute.\n");
          LADA_ASSERT(vec_node->Attribute("y"), "Missing y attribute.\n");
          LADA_ASSERT(vec_node->Attribute("z"), "Missing z attribute.\n");
          MLCluster::Spin const spin = 
          {
            0,
            math::rVector3d
            (
              boost::lexical_cast<types::t_real>(vec_node->Attribute("x")), 
              boost::lexical_cast<types::t_real>(vec_node->Attribute("y")), 
              boost::lexical_cast<types::t_real>(vec_node->Attribute("z"))  
            ) 
          };
          
          // avoid origin.
          if( math::is_null( spin.pos.squaredNorm() ) ) continue;
          cluster.push_back(spin);
        }

        clusters.init(cluster, clusters.eci);
        result->push_back(clusters);
      }
      return result;
    }  

    boost::shared_ptr<t_MLClusterClasses> load_new(TiXmlElement const &_node)
    {
      
      TiXmlElement const * clusters_node = opt::find_node(_node, "Clusters");
      LADA_DOASSERT( clusters_node, "Could not find Clusters node in input.\n" );
      TiXmlElement const * cluster_node = clusters_node->FirstChildElement("Cluster");
      LADA_DOASSERT( cluster_node, "Could not find Clusters node in input.\n" );
       
      boost::shared_ptr<t_MLClusterClasses> result(new t_MLClusterClasses);
      for(; cluster_node; cluster_node = cluster_node->NextSiblingElement("Cluster") )
      {
        result->push_back(MLClusters());
        result->back().load( *cluster_node );
      }
      return result;
    } 
  } // namespace CE
}
