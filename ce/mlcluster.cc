//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <boost/filesystem/operations.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/regex.hpp>

#include <crystal/neighbors.h>
#include <crystal/symmetry_operator.h>
#include <crystal/which_site.h>
#include <opt/types.h>
#include <opt/fuzzy.h>
#include <opt/atat.h>
#include <atat/findsym.h>
#include <atat/xtalutil.h>

#include "cluster.h"


namespace LaDa
{
  namespace CE 
  {

    void MLCluster :: apply_symmetry( Crystal::SymmetryOperator const &_op )
    {
      if( is_null ) return;
      // finds new origin.
      atat::rVector3d const transformed_pos(_op(lattice.sites[origin].pos)); 
      LADA_ASSERT( which_site(transformed_pos, inv_cell, lattice.sites) != -1, 
                   "Rotation is not of matrix.\n" );
      origin = std:abs(which_site(transformed_pos, inv_cell, lattice.sites));

      if ( vectors.size() < 1 ) return;

      // finds new vectors and sites.
      LADA_ASSERT(Crystal::Structure::lattice != NULL, "No lattice is set.\n");
      Crystal::Lattice const &lattice = *Crystal::Structure::lattice;

      LADA_ASSERT(origin < lattice.sites.size(), "Index out of range.\n");
      LADA_ASSERT(vectors.size() == sites.size(), "Incoherent number of sites and vectors.\n");

      atat::rMatrix3d const inv_cell( !lattice.cell );

      t_Positions :: iterator i_vec = vectors.begin();
      t_Positions :: iterator const i_last = vectors.end();
      t_Sites :: iterator i_site = sites.begin();
      for(; i_vec != i_last; ++i_vec, ++i_site)
      {
        *i_vec = _op(*i_vec);
        *i_site = which_site( *i_vec, inv_cell, lattice.sites );
      }

      // shifts back to origin if necessary.
      if( not Fuzzy::is_zero(atat::norm2(transformed_pos)) ) 
        for(i_vec = vectors.begin(); i_vec != i_last; ++i_vec ) *i_vec -= transformed_pos; 
    }


    std::ostream& operator<<( std::ostream &_stream,  MLCluster const &_cls )
    {
      _stream << std::fixed << std::setprecision(5) << std::setw(9);
      _stream << "Cluster@" << _clsorigin;
      
      if (_cls.is_null) return _stream << " J0\n";
      if (_cls.vectors.size() == 0 ) return << " J1\n";

      LADA_ASSERT(_cls.vectors.size() == _cls.sites.size(), "Incoherent sites and vectors.\n");
      t_Sites :: const_iterator i_site = _cls.sites.begin();
      t_Positions :: const_iterator i_vec = _cls.vectors.begin();
      t_Positions :: const_iterator const i_last = _cls.vectors.end();
      for ( ; i_vec != i_last; ++i_vec, ++i_site)
        _stream << " @" << *i_site
                << " " << ( *i_vec )(0) 
                << " " << ( *i_vec )(1) 
                << " " << ( *i_vec )(2) 
                << "\n";
      return _stream;
    }

    bool Cluster :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement *child;
      types::t_real d; atat::rVector3d vec;

      vectors.clear();
      sites.clear();

      is_null = not _node.Attribute("site"); 
      if( is_null )  return true;

      origin = boost::lexical_cast<size_t>( child->Attribute("site") );
      child = _node.FirstChildElement( "spin" );
      for ( ; child; child=child->NextSiblingElement( "spin" ) )
      {
        LADA_ASSERT( child->Attribute("x"), "Could not find x attribute.\n" );
        LADA_ASSERT( child->Attribute("y"), "Could not find y attribute.\n" );
        LADA_ASSERT( child->Attribute("z"), "Could not find z attribute.\n" );
        LADA_ASSERT( child->Attribute("site"), "Could not find site attribute.\n" );
        vec(0) = boost::lexical_cast<types::t_real>(child->Attribute("x"));
        vec(1) = boost::lexical_cast<types::t_real>(child->Attribute("y"));
        vec(2) = boost::lexical_cast<types::t_real>(child->Attribute("z"));
        sites.push_back( boost::lexical_cast<size_t>( child->Attribute("site") ) );
        vectors.push_back(vec);
      }

      return true;
    }

    // returns true if second position in cluster is equal to pos.
    struct CompPairs
    {
      atat::rVector3d const pos;
      CompPairs( atat::rVector3d const &_a ) : pos(_a) {}
      CompPairs( CompPairs const &_a ) : pos(_a.pos) {}
      bool operator()(Cluster const& _a) const
        { return Fuzzy::is_zero(atat::norm2(_a.vectors[1]-pos)); }
    };

   
    void read_clusters( const Crystal::Lattice &_lat, 
                        const boost::filesystem::path &_path, 
                        std::vector< std::vector< Cluster > > &_out,
                        const std::string &_genes )
    {
      __DEBUGTRYBEGIN

      namespace fs = boost::filesystem;  
      __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                  _path << " is neither a regulare file nor a system link.\n" )
      std::ifstream file( _path.string().c_str(), std::ifstream::in );
      std::string line;

      Cluster cluster;
      size_t i(0);
      while( read_cluster( _lat, file, cluster ) )
      {
        // This is a hack to avoid inputing Zhe's tetragonal pair terms.
        if( cluster.size() == 2 ) continue;
        if( not _genes.empty() )
        {
          __DOASSERT( i >= _genes.size(), "Gene size and jtypes are inconsistent.\n" )
          if( _genes[i] == '0' ) { ++i; continue; }
        }
        ++i;
        // creates a new class from cluster.
        _out.push_back( std::vector< Cluster >( 1, cluster ) );
        // adds symmetrically equivalent clusters.
        add_equivalent_clusters( _lat, _out.back() );
      }
      if( not _genes.empty() )
        __DOASSERT( i != _genes.size(), "Gene size and jtypes are inconsistent.\n" )

      __DEBUGTRYEND(,"Could not read clusters from input.\n" )
    }
    bool read_cluster( const Crystal::Lattice &_lat, 
                       std::istream & _sstr,
                       Cluster &_out )
    {
      __DEBUGTRYBEGIN
      std::string line;
      // bypass comments until a name is found.
      const boost::regex comment("^(\\s+)?#");
      boost::match_results<std::string::const_iterator> what;

      do
      { // check for comment.
        do { std::getline( _sstr, line ); }
        while(     boost::regex_search( line, what, comment )
               and ( not _sstr.eof() ) );

        if( _sstr.eof() ) return false;
      } while( not (    line.find( 'B' ) != std::string::npos 
                     or line.find( 'J' ) != std::string::npos  ) );

      // now read cluster.
      // second line should contain order and multiplicity.
      { // check for comment.
        std::getline( _sstr, line );
        while(     boost::regex_search( line, what, comment )
               and ( not _sstr.eof() ) )
          std::getline( _sstr, line ); 
        if( _sstr.eof() ) return false;
      }
      types::t_unsigned order;
      types::t_unsigned multiplicity;
      { // should contain the order and the multiplicity.
        std::istringstream sstr( line ); 
        sstr >> order >> multiplicity;
        __DOASSERT( sstr.bad(), "Error while reading figure.\n" );
      }

      // creates cluster.
      Cluster cluster;
      for(; order; --order)
      {
        { // check for comment.
          std::getline( _sstr, line ); 
          while(     boost::regex_search( line, what, comment )
                 and ( not _sstr.eof() ) )
            std::getline( _sstr, line ); 
          __DOASSERT( _sstr.eof(), "Unexpected end of file.\n" )
        }
        __DOASSERT( _sstr.eof(),    "Unexpected end-of-file.\n" )
        std::istringstream sstr(line);
        types::t_real x, y, z;
        sstr >> x >> y >> z;
        cluster.vectors.push_back( atat::rVector3d( x * 5e-1, y * 5e-1, z * 5e-1 ) );
      }
      _out = cluster;
      return true;
      __DEBUGTRYEND(,"Could not read clusters from input.\n" )
    }

    void  add_equivalent_clusters( const Crystal::Lattice &_lat, 
                                   std::vector< Cluster > &_out )
    {
      typedef std::vector< Cluster > t_Clusters;
      atat::rMatrix3d const &cell    = _lat.cell;
      atat::rMatrix3d const inv_cell = !cell;
      // new clusters will be added directly to _out
      t_Clusters old_cluster_list = _out;
    
      t_Clusters :: iterator i_old_list_last = old_cluster_list.end();
      t_Clusters :: iterator i_old_list = old_cluster_list.begin();
      for (; i_old_list != i_old_list_last ; ++i_old_list )
      {
        typedef std::vector<Crystal::SymmetryOperator> :: const_iterator t_cit;
        t_cit i_op = _lat.space_group.begin();
        t_cit const i_op_end = _lat.space_group.end();
        for (; i_op != i_op_end; ++i_op )
        {
          if( not Fuzzy::is_zero(atat::norm2(i_op->trans)) ) continue;
          // initialize a new cluster to object pointed by i_old_list
          Cluster transfo_cluster( *i_old_list );
          
          // transforms cluster according to symmetry group
          transfo_cluster.apply_symmetry(*i_op);
    
          // checks wether transformed cluster is in Lamarck::clusters
          t_Clusters :: iterator i_cluster  = _out.begin();
          t_Clusters :: iterator i_last = _out.end();
          for ( ; i_cluster != i_last ; ++i_cluster)
            if ( transfo_cluster.equivalent_mod_cell(*i_cluster, inv_cell) ) 
              break;
    
          // if it isn't, adds cluster to clusters
          if (i_cluster != i_last ) continue;
          _out.push_back( transfo_cluster );
        }
      }
    }
    
    // == for atat::rVector3d
    struct cmp
    {
      cmp(atat::rVector3d const &_c) : pos(_c) {}
      cmp(cmp const &_c): pos(_c.pos) {}
      bool operator()(atat::rVector3d const &_a) const
      {
        if( not Fuzzy::is_zero( pos[0] - _a[0] ) ) return false;
        if( not Fuzzy::is_zero( pos[1] - _a[1] ) ) return false;
        return  Fuzzy::is_zero( pos[2] - _a[2] );
      }
      atat::rVector3d const &pos;
    };

    bool Cluster::operator==( Cluster const & _c ) const
    {
      if( vectors.size() != _c.vectors.size() ) return false;

      foreach( atat::rVector3d const &pos, _c.vectors )
        if( vectors.end() == std::find_if(vectors.begin(), vectors.end(), cmp(pos)) ) 
          return false;
      return true;
    }
  } // namespace CE
}
