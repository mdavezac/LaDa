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

#include <crystal/symmetry_operator.h>
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

    void Cluster :: apply_symmetry( Crystal::SymmetryOperator const &_op )
    {
      if ( vectors.size() < 1 ) return;
      std::vector<atat::rVector3d> :: iterator i_vec = vectors.begin();
      std::vector<atat::rVector3d> :: iterator const i_last = vectors.end();
      for(; i_vec != i_last; ++i_vec) *i_vec = _op(*i_vec);

      if( Crystal::Structure::lattice == NULL ) return;

      i_vec = vectors.begin();
      atat::rVector3d shift( (!Crystal::Structure::lattice->cell) * (*i_vec) );
      shift(0) -= std::floor(shift(0)); if( Fuzzy::is_zero(shift(0)-1e0) ) shift(0) = 0e0;
      shift(1) -= std::floor(shift(1)); if( Fuzzy::is_zero(shift(1)-1e0) ) shift(1) = 0e0;
      shift(2) -= std::floor(shift(2)); if( Fuzzy::is_zero(shift(2)-1e0) ) shift(2) = 0e0;
      shift = Crystal::Structure::lattice->cell * shift - (*i_vec);
      if( Fuzzy::is_zero(atat::norm2(frac)) ) return;
      for(; i_vec != i_vec_end; ++i_vec ) *i_vec -= frac; 
    }

    bool Cluster :: equivalent_mod_cell( Cluster &_cluster, const atat::rMatrix3d &_icell) 
    {
      namespace bl = boost::lambda;
      if ( vectors.size() != _cluster.vectors.size() ) return false;
      if ( vectors.size() == 0 ) return true;

      // Gets fractional coordinates.
      std::vector<atat::rVector3d> fractionals;
      std::vector<atat::rVector3d> :: iterator i_vec = _clusters.vectors.begin();
      std::vector<atat::rVector3d> :: iterator i_vec_last = _clusters.vectors.end();
      for (; i_vec != i_vec_last; ++i_vec)
      {
        atat::rVector3d frac(  _icell * (*i_vec) );
        frac(0) -= std::floor(frac(0)); if( Fuzzy::is_zero(frac(0)-1e0) ) frac(0) = 0e0;
        frac(1) -= std::floor(frac(1)); if( Fuzzy::is_zero(frac(1)-1e0) ) frac(1) = 0e0;
        frac(2) -= std::floor(frac(2)); if( Fuzzy::is_zero(frac(2)-1e0) ) frac(2) = 0e0;
        fractionals.push_back(frac);
      }
     

      i_vec = vectors.begin();
      i_vec_end = vectors.end();
      for (; i_vec != i_vec_last; ++i_vec)
      {
        atat::rVector3d frac(  _icell * (*i_vec) );
        frac(0) -= std::floor(frac(0)); if( Fuzzy::is_zero(frac(0)-1e0) ) frac(0) = 0e0;
        frac(1) -= std::floor(frac(1)); if( Fuzzy::is_zero(frac(1)-1e0) ) frac(1) = 0e0;
        frac(2) -= std::floor(frac(2)); if( Fuzzy::is_zero(frac(2)-1e0) ) frac(2) = 0e0;
        if( fractionals.end() == std::find(fractionals.begin(), fractionals.end(), frac) )
          return false;
      }
      return true;
    }

    void Cluster :: print_out (  std::ostream &stream)  const
    {
      stream << std::fixed << std::setprecision(5) << std::setw(9);
      stream << " Cluster: " << eci << std::endl;
      
      if (vectors.size() == 0 )
      {
        stream << " No position" << std::endl;
        return;
      }

      std::vector<atat::rVector3d> :: const_iterator i_vec = vectors.begin();
      std::vector<atat::rVector3d> :: const_iterator i_last = vectors.end();
      
      for ( ; i_vec != i_last; ++i_vec)
        stream << " " << ( *i_vec )(0) 
               << " " << ( *i_vec )(1) 
               << " " << ( *i_vec )(2) 
               << std::endl;
    }

    bool Cluster :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement *child;
      types::t_real d; atat::rVector3d vec;

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


    void create_pairs( const Crystal :: Lattice &_lat,
                       types::t_unsigned _max_neigh,
                       std::vector< std::vector<Cluster> > &_out )
    {
      __DEBUGTRYBEGIN
      typedef std::vector< std::vector< Cluster > > t_EquivClusters;
      if( _max_neigh == 0 ) return;
      
      namespace bl = boost::lambda;
      // Creates a vector of all positions within a cube \a _max_neigh big.
      std::vector< atat::rVector3d > nn;
      nn.reserve( _max_neigh*_max_neigh*_max_neigh );
      for( types::t_real x(-types::t_real(_max_neigh));
           Fuzzy::leq(x, (types::t_real)_max_neigh); x+=1e0 )
        for( types::t_real y(-types::t_real(_max_neigh));
             Fuzzy::leq(y, (types::t_real)_max_neigh); y+=1e0 )
          for( types::t_real z(-types::t_real(_max_neigh));
               Fuzzy::leq(z, (types::t_real)_max_neigh); z+=1e0 )
            if(    Fuzzy::neq( x, types::t_real(0) )
                or Fuzzy::neq( y, types::t_real(0) )
                or Fuzzy::neq( z, types::t_real(0) ) )
            nn.push_back( _lat.cell * atat::rVector3d( x, y, z ) );

      // Sorts these positions by size.
      types::t_real (*ptr_norm2)( const atat::FixedVector<types::t_real,3>& ) 
                           = &atat::norm2<types::t_real, 3>;
      std::sort
      ( 
        nn.begin(), nn.end(), 
        bl::bind
        (
          std::less<types::t_real>(),
          bl::bind<types::t_real>( std::ptr_fun( ptr_norm2 ), bl::_1 ),
          bl::bind<types::t_real>( std::ptr_fun( ptr_norm2 ), bl::_2 ) 
        )
      );

      // finds _max_neigh neighbor.
      std::vector< atat :: rVector3d > :: iterator n_it = nn.begin();
      __DODEBUGCODE( std::vector< atat :: rVector3d > :: iterator n_it_end = nn.end(); )
      types::t_real norm = atat::norm2( *n_it ); 
      for( types::t_unsigned N( _max_neigh ); N != 0; ++n_it )
      {
        __ASSERT( n_it == n_it_end, "Iterator out of range with N=" << N << ".\n" )
        if( Fuzzy::eq( norm, atat::norm2( *n_it ) ) ) continue;
        --N;
        norm = atat::norm2( *n_it );
      }
      __ASSERT( n_it - nn.begin() -1 <= 0, "Requesting vector of size <= 0.\n" )
      nn.resize( n_it - nn.begin() - 1 );
      
      // Create classes of equivalent positions.
      std::vector< std::vector< atat::rVector3d > > classes, onesize;
      //    -- first adds first (null) position.
      onesize.push_back( std::vector< atat::rVector3d >( 1, nn.front() ) );
      n_it = nn.begin() + 1;
      __NDEBUGCODE(std::vector< atat :: rVector3d > :: iterator)
        n_it_end = nn.end();
      norm = atat::norm2( *n_it );
      //    -- then goes across all neighboring atoms and organizes in classes.
      for(; n_it != n_it_end; ++n_it )
      {
        if( Fuzzy::neq( norm, atat::norm2( *n_it ) ) )
        {
          std::copy( onesize.begin(), onesize.end(),
                     std::back_inserter( classes ) );
          onesize.clear();
          onesize.push_back( std::vector< atat::rVector3d >( 1, *n_it ) );
          norm = atat::norm2( *n_it );
          continue;
        }
        // looks for position in class i_found.
        std::vector< std::vector< atat::rVector3d > > :: iterator
          i_found = onesize.begin();
        std::vector< std::vector< atat::rVector3d > > :: iterator
          i_found_end = onesize.end();
        for(; i_found != i_found_end; ++i_found )
        {
          // checks if position is equivalent to something in class i_found.
          std::vector< atat::rVector3d > :: const_iterator i_equiv = i_found->begin();
          std::vector< atat::rVector3d > :: const_iterator i_notequiv = i_found->end();
          for(; i_equiv != i_notequiv; ++i_equiv )
            if( _lat.equiv_by_point_group( *i_equiv, *n_it ) ) break;
          // there are no equivalents, continue;
          if( i_equiv == i_notequiv ) continue;
          // adds to i_found and breaks.
          i_found->push_back( *n_it );
          break;
        }
        if( i_found != onesize.end() ) continue;
        // Current position is new. Adds a class to onesize.
        onesize.push_back( std::vector< atat::rVector3d >( 1, *n_it ) );
      }
      if( onesize.size() ) std::copy( onesize.begin(), onesize.end(),
                                      std::back_inserter( classes ) );
      onesize.clear();

      Cluster cluster;
      cluster.vectors.resize(2);
      cluster.vectors[0] = atat::rVector3d(0,0,0); 
      cluster.eci = 0e0;
      std::vector< std::vector< atat::rVector3d > > :: const_iterator
        i_class = classes.begin();
      std::vector< std::vector< atat::rVector3d > > :: const_iterator
        i_class_end = classes.end();
      for(; i_class != i_class_end; ++i_class )
      {
        _out.resize( _out.size() + 1 );
        t_EquivClusters :: value_type& back = _out.back();
        std::vector< atat::rVector3d > :: const_iterator i_pos = i_class->begin();
        std::vector< atat::rVector3d > :: const_iterator i_pos_end = i_class->end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          cluster.vectors.back() = *i_pos;
          back.push_back( cluster );
        }
      }

      __DEBUGTRYEND(,"Could not create pair clusters.\n" )
    }
   
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
        details::add_equivalent_clusters( _lat, _out.back() );
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
      LADA_ASSERT( _site < _lat.sites.size(), "Index out of range.\n" );
      typedef std::vector< Cluster > t_Clusters;
      atat::rMatrix3d const &cell    = _lat.cell;
      atat::rMatrix3d const inv_cell = !cell;
      // new clusters will be added directly to _out
      t_Clusters old_cluster_list = _out;
    
      t_Clusters :: iterator i_old_list_last = old_cluster_list.end();
      t_Clusters :: iterator i_old_list = old_cluster_list.begin();
      for (; i_old_list != i_old_list_last ; ++i_old_list )
      {
        std::vector<SymmetryOperator> :: const_iterator i_op = _lat.space_group.begin();
        std::vector<SymmetryOperator> :: const_iterator const i_op_end = _lat.space_group.end();
        for (; i_op != i_op_end; ++i_op )
        {
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
    
  } // namespace CE
}
