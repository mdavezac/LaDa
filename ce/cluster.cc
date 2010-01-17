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
#include <opt/types.h>
#include <math/fuzzy.h>
#include <math/compare_norms.h>
#include <math/misc.h>

#include "cluster.h"


namespace LaDa
{
  namespace CE 
  {

    void Cluster :: apply_symmetry( Crystal::SymmetryOperator const &_op )
    {
      if ( vectors.size() < 1 ) return;
      std::vector<math::rVector3d> :: iterator i_vec = vectors.begin();
      std::vector<math::rVector3d> :: iterator const i_last = vectors.end();
      for(; i_vec != i_last; ++i_vec) *i_vec = _op(*i_vec);

      if( Crystal::Structure::lattice == NULL ) return;

      i_vec = vectors.begin();
      math::rVector3d shift( (!Crystal::Structure::lattice->cell) * (*i_vec) );
      shift(0) -= std::floor(shift(0)); if( math::is_zero(shift(0)-1e0) ) shift(0) = 0e0;
      shift(1) -= std::floor(shift(1)); if( math::is_zero(shift(1)-1e0) ) shift(1) = 0e0;
      shift(2) -= std::floor(shift(2)); if( math::is_zero(shift(2)-1e0) ) shift(2) = 0e0;
      shift = Crystal::Structure::lattice->cell * shift - (*i_vec);
      if( math::is_zero(shift.squaredNorm()) ) return;
      for(; i_vec != i_last; ++i_vec ) *i_vec -= shift; 
    }

    bool Cluster :: are_periodic_images( Cluster &_cluster, const math::rMatrix3d &_icell) 
    {
      if ( vectors.size() != _cluster.vectors.size() ) return false;
      if ( vectors.size() == 0 ) return true;
      
      std::vector<math::rVector3d> :: iterator i_vec = vectors.begin();
      std::vector<math::rVector3d> :: iterator i_vec_last = vectors.end();
      for (; i_vec != i_vec_last; ++i_vec)
      {
        if ( not math::are_periodic_images( *_cluster.vectors.begin(), *i_vec, _icell) ) continue;
      
        math::rVector3d shift = (*i_vec) - *(_cluster.vectors.begin());
        std::vector<math::rVector3d> :: iterator is_found;
        
        // search for _cluster  vector such that|| (*i_vec-shift) - *i_equiv ||  > zero_tolerance
        std::vector<math::rVector3d> :: iterator i2_vec = vectors.begin();
        for(; i2_vec != i_vec_last; ++i2_vec)
        {
          is_found  = std::find_if( _cluster.vectors.begin(),
                                    _cluster.vectors.end(),
                                    math::CompareNorms((*i2_vec - shift).eval()) );
      
          if ( is_found == _cluster.vectors.end() ) break;
        }
                                      
        // if all match, return true
        if ( i2_vec == i_vec_last ) return true; 
        
      }
      return false;
    }

    void Cluster :: print_out (  std::ostream &stream)  const
    {
      stream << std::fixed << std::setprecision(5) << std::setw(9);
      stream << "Cluster: " << eci << "\n";
      
      if (vectors.size() == 0 )
      {
        stream << " No position" << std::endl;
        return;
      }

      std::vector<math::rVector3d> :: const_iterator i_vec = vectors.begin();
      std::vector<math::rVector3d> :: const_iterator i_last = vectors.end();
      
      for ( ; i_vec != i_last; ++i_vec)
        stream << " " << ( *i_vec )(0) 
               << " " << ( *i_vec )(1) 
               << " " << ( *i_vec )(2) 
               << "\n";
    }

    bool Cluster :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement *child;
      types::t_real d; math::rVector3d vec;

      _node.Attribute("eci", &eci);
      vectors.clear();
      child = _node.FirstChildElement( "spin" );
      for ( ; child; child=child->NextSiblingElement( "spin" ) )
      {
        d = 1.0 ; child->Attribute("x", &d); vec(0) = d;
        d = 1.0 ; child->Attribute("y", &d); vec(1) = d;
        d = 1.0 ; child->Attribute("z", &d); vec(2) = d;
        vectors.push_back(vec);
      }

      return true;
    }

    // returns true if second position in cluster is equal to pos.
    struct CompPairs
    {
      math::rVector3d const pos;
      CompPairs( math::rVector3d const &_a ) : pos(_a) {}
      CompPairs( CompPairs const &_a ) : pos(_a.pos) {}
      bool operator()(Cluster const& _a) const
        { return math::is_zero( (_a.vectors[1]-pos).squaredNorm() ); }
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
        cluster.vectors.push_back( math::rVector3d( x * 5e-1, y * 5e-1, z * 5e-1 ) );
      }
      _out = cluster;
      return true;
      __DEBUGTRYEND(,"Could not read clusters from input.\n" )
    }

    void  add_equivalent_clusters( const Crystal::Lattice &_lat, 
                                   std::vector< Cluster > &_out )
    {
      typedef std::vector< Cluster > t_Clusters;
      math::rMatrix3d const &cell    = _lat.cell;
      math::rMatrix3d const inv_cell = cell.inverse();
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
          if( not math::is_zero(i_op->trans.squaredNorm()) ) continue;
          // initialize a new cluster to object pointed by i_old_list
          Cluster transfo_cluster( *i_old_list );
          
          // transforms cluster according to symmetry group
          transfo_cluster.apply_symmetry(*i_op);
    
          // checks wether transformed cluster is in Lamarck::clusters
          t_Clusters :: iterator i_cluster  = _out.begin();
          t_Clusters :: iterator i_last = _out.end();
          for ( ; i_cluster != i_last ; ++i_cluster)
            if ( transfo_cluster.are_periodic_images(*i_cluster, inv_cell) ) 
              break;
    
          // if it isn't, adds cluster to clusters
          if (i_cluster != i_last ) continue;
          _out.push_back( transfo_cluster );
        }
      }
    }
    
    // == for math::rVector3d
    struct cmp
    {
      cmp(math::rVector3d const &_c) : pos(_c) {}
      cmp(cmp const &_c): pos(_c.pos) {}
      bool operator()(math::rVector3d const &_a) const
      {
        if( not math::is_zero( pos[0] - _a[0] ) ) return false;
        if( not math::is_zero( pos[1] - _a[1] ) ) return false;
        return  math::is_zero( pos[2] - _a[2] );
      }
      math::rVector3d const &pos;
    };

    bool Cluster::operator==( Cluster const & _c ) const
    {
      if( vectors.size() != _c.vectors.size() ) return false;

      foreach( math::rVector3d const &pos, _c.vectors )
        if( vectors.end() == std::find_if(vectors.begin(), vectors.end(), cmp(pos)) ) 
          return false;
      return true;
    }
  } // namespace CE
}
