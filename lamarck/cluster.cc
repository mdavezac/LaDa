//
//  Version: $Id$
//
#include <string>
#include <iomanip>
#include <algorithm>
#include <boost/filesystem/operations.hpp>
#include <boost/lambda/bind.hpp>

#include <opt/types.h>
#include <opt/fuzzy.h>
#include <opt/atat.h>
#include <atat/findsym.h>
#include <atat/xtalutil.h>

#include "cluster.h"


namespace CE {

  void Cluster :: apply_symmetry( const atat::rMatrix3d &_op,
                                  const atat::rVector3d &_trans    ) 
  {
    if ( vectors.size() < 1 ) return;
    std::vector<atat::rVector3d> :: iterator i_vec = vectors.begin();
    std::vector<atat::rVector3d> :: iterator i_last = vectors.end();
    for(; i_vec != i_last; ++i_vec)
      *i_vec = _op * (*i_vec) + _trans;

    i_vec = vectors.begin();
    atat::rVector3d translate = *i_vec;
    for(; i_vec != i_last; ++i_vec)
      *i_vec -= translate;
  }

  bool Cluster :: equivalent_mod_cell( Cluster &_cluster, const atat::rMatrix3d &_icell) 
  {
    if ( vectors.size() != _cluster.vectors.size() ) return false;
    if ( vectors.size() == 0 ) return true;

    std::vector<atat::rVector3d> :: iterator i_vec = vectors.begin();
    std::vector<atat::rVector3d> :: iterator i_vec_last = vectors.end();
    for (; i_vec != i_vec_last; ++i_vec)
    {
      if ( not atat::equivalent_mod_cell( *_cluster.vectors.begin(), *i_vec, _icell) ) continue;

      atat::rVector3d shift = (*i_vec) - *(_cluster.vectors.begin());
      std::vector<atat::rVector3d> :: iterator is_found;
      
      // search for _cluster  vector such that|| (*i_vec-shift) - *i_equiv ||  > zero_tolerance
      std::vector<atat::rVector3d> :: iterator i2_vec = vectors.begin();
      for(; i2_vec != i_vec_last; ++i2_vec)
      {
        is_found  = std::find_if( _cluster.vectors.begin(),
                                  _cluster.vectors.end(),
                                  atat::norm_compare( *i2_vec - shift ) );

        if ( is_found == _cluster.vectors.end() ) break;
      }
                                    
      // if all match, return true
      if ( i2_vec == i_vec_last ) return true; 
      
    }
    return false;
  }

  void Cluster :: print_out (  std::ostream &stream)  const
  {
    stream << std::fixed << std::setprecision(2) << std::setw(6);
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


  void find_all_clusters( const Crystal :: Lattice &_lat,
                          types::t_unsigned _max_neigh,
                          types::t_unsigned _maxN,
                          std::vector< std::vector<Cluster> > &_out )
  {
    __DEBUGTRYBEGIN

    namespace bl = boost::lambda;
    // Creates a vector of all positions within a cube \a _max_neigh big.
    std::vector< atat::rVector3d > nn;
    nn.reserve( _max_neigh**3 );
    for( types::t_real x(-_max_neigh);
         Fuzzy::lt(x, (types::t_real)_max_neigh); x+=1e0 )
      for( types::t_real y(-_max_neigh);
           Fuzzy::lt(y, (types::t_real)_max_neigh); y+=1e0 )
        for( types::t_real z(-_max_neigh);
             Fuzzy::lt(z, (types::t_real)_max_neigh); z+=1e0 )
          nn.push_back( _lat.cell * atat::rVector3d( x, y, z ) );

    // Sorts these positions by size.
    std::sort
    ( 
      nn.begin(), nn.end(), 
      bl::bind( &atat::norm2, bl::_1 ) <= bl::bind( &atat::norm2, bl::_2 )
    );

    // finds _max_neigh neighbor.
    std::vector< atat :: rVector3d > :: iterator n_it = nn.begin();
    __DODEBUGCODE( std::vector< atat :: rVector3d > :: iterator n_it_end = nn.end(); )
    types::t_real norm = atat::norm2( *n_it );
    for( types::t_unsigned N( _max_neigh ); N != 0; ++n_it )
    {
      __ASSERT( n_it != n_it_end, "Iterator out of range.\n" )
      if( Fuzzy::eq( norm, atat::norm2( *n_it ) ) ) continue;
      --N;
      norm = atat::norm2( *n_it );
    }
    __ASSERT( n_it - nn.begin() -1 <= 0, "Requesting vector of size <= 0.\n" )
    nn.resize( n_it - nn.begin() - 1 );
    
    // Create classes of equivalent positions.
    std::vector< vector< atat::rVector3d > > classes, onesize;
    //    -- first adds first (null) position.
    onesize.push_back( std::vector< atat::rVector3d >( 1, nn.front() ) );
    n_it = nn.begin(); ++n_it;
    __NDEBUGCODE(std::vector< atat :: rVector3d > :: iterator)
      n_it_end = nn.end();
    norm = 0e0;
    //    -- then goes across all neighboring atoms and organizes in classes.
    for(; n_it != n_it_end; ++n_it )
    {
      if( Fuzzy::neq( norm, atat::norm2( *n_it ) ) )
      {
        classes.append( onesize );
        onesize.clear();
        onesize.push_back( std::vector< atat::rVector3d >( 1, *n_it ) );
        norm = atat::norm2( *n_it );
        continue;
      }
      // looks for position in class i_found.
      std::vector< vector< atat::rVector3d > > :: iterator i_found = onesize.begin();
      std::vector< vector< atat::rVector3d > > :: iterator i_found_end = onesize.end();
      for(; i_found != i_found_end; ++i_found )
      {
        std::vector< atat::rVector3d > i_equiv;
        // checks if position is equivalent to something in class i_found.
        i_equiv = std::find_if
                  (
                    i_found->begin(), i_found->end(),
                    bl::bind( &Lattice::pos_are_equiv, &lattice, bl::_1, bl::_2 ) 
                  ); 
        // there are no equivalents, continue;
        if( i_equiv == i_found->end() ) continue;
        // adds to i_found and breaks.
        i_found->add( *n_it );
        break;
      }
      if( i_found != onesize.end() ) continue;
      // Current position is new. Adds a class to onesize.
      onesize.push_back( std::vector< atat::rVector3d >( 1, *n_it ) );
    }
    if( onesize.size() ) classes.append( onesize );
    onesize.clear();

    // Creates classes of clusters, by taking one of each equivalent position class.
    opt::Ndim_Iterator<size_t, std::less<size_t> > ndim;
    for(size_t N( _maxN -1 ); N; --N )
      ndim.add( 0, classes.size() );

    std::vector< std::vector< Cluster > > newclusters;
    do
    {
      // makes sure indices are ordered k > .. > j > i.
      bool notordered = false;
      do
      {
        for(size_t N( _maxN - 1 ); N-1; --N )
        {
          if( ndim.acess[N] > ndim.access[N-1] ) continue;
          notordered = true;
          break;
        }
      }
      while( notordered );

      // now creates two body figures.
      Cluster cluster;
      cluster.vectors.resize(2, atat::rVector3d(0,0,0) );
      std::vector< std::vector< Cluster > > new_clusters(1);
      std::vector< atat::rVector3d > :: const_iterator i_pos = classes[ ndim.access(0) ].begin();
      std::vector< atat::rVector3d > :: const_iterator i_pos_end = classes[ ndim.access(0) ].end();
      for(; i_pos != i_pos_end; ++i_pos )
      {
        cluster.vectors.back() = *i_pos;
        new_clusters.front().push_back( cluster ); 
      }

      // then adds higher order clusters, if required.
      if( _maxN == 2 ) continue;
      for( size_t i(1); i < _maxN; ++i )
      {
        std::vector< Cluster > :: const_iterator i_cluster = new_clusters.back().begin();
        std::vector< Cluster > :: const_iterator i_cluster_end = new_clusters.back().end();
        new_clusters.resize( new_clusters.size() + 1 );
        for(; i_cluster != i_cluster_end; ++i_cluster )
        {
          cluster = *i_cluster;
          cluster.vectors.push_back( atat::rVector3d(0,0,0) );
          i_pos = classes[ ndim.access( i ) ].begin();
          i_pos_end = classes[ ndim.access( i ) ].end();
          for(; i_pos != i_pos_end; ++i_pos )
          {
            cluster.vectors.back() = *i_pos;
            new_clusters.back().push_back( cluster );
          }
        }
      } // end of higher order clusters.

    } while( ++ndim );

    // Finally, appends new cluster classes to output vector.
    // At this point, we could verify that the classes do not already exist in
    // the output.
    _out.append( new_clusters );

    __DEBUGTRYEND(,"Could not create clusters.\n" )
  }
 
  void read_clusters( Crystal::Lattice &_lat, 
                      boost::filesystem::path &_path, 
                      std::vector< std::vector< Cluster > > &_out );
  {
    __DEBUGTRYBEGIN

    namespace fs = boost::filesystem;  
    __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
    __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                _path << " is neither a regulare file nor a system link.\n" )
    std::ifstream file( fullpath.string().c_str(), std::ifstream::in );
    std::string line;
    size_t i(0);

    while( std::getline( file, line ) )
    {
      ++i;
      // first line should be a name. Ignore it.
      std::getline( file, line ); ++i;
      __DOASSERT( not file.eof(),    "Unexpected end of file at line "
                                  << i << " in " << _path << ".\n" )
      std::istringstream sstr( line ); // should contain the order and the multiplicity.
      types::t_unsigned order;
      types::t_unsigned multiplicity;
      sstr >> order >> multiplicity;
      __DOASSERT( sstr.bad(), "Error while reading figure.\n" );
      // creates cluster.
      Cluster cluster;
      for(; order; --order)
      {
        std::getline( file, line ); ++i;
        __DOASSERT( not file.eof(),    "Unexpected end of file at line "
                                    << i << " in " << _path << ".\n" )
        sstr.str(line);
        atat::rVector3d pos;
        sstr >> pos[0] >> pos[1] >> pos[2];
        cluster.vectors.push_back( pos );
      }
      // creates a new class from cluster.
      _out.push_back( std::vector< Cluster >( 1, cluster ) );
      // adds symmetrically equivalent clusters.
      details::add_equivalent_clusters( _lat, _out.back() );
    }

    __DEBUGTRYEND(,"Could not read clusters from input.\n" )
  }

  namespace details
  {
    // finds all clusters, including symmetric equivalents
    // starting from cluster included in Lamarck::clusters
    // results are stored in a new Lamarck::clusters
    void  add_equivalent_clusters( Crystal::Lattice &_lat, 
                                   std::vector< Cluster > &_out  ) 
    {
      typedef std::vector< Cluster > t_Clusters;
      const atat::rMatrix3d &cell            = _lat.cell;
      const atat::rMatrix3d inv_cell         = !cell;
      const atat::Array<atat::rMatrix3d> &point_op = _lat.space_group.point_op;
      const atat::Array<atat::rVector3d> &trans    = _lat.space_group.trans;
      // new clusters will be added directly to _out
      t_Clusters old_cluster_list = _out;
  
      t_Clusters :: iterator i_old_list_last = old_cluster_list.end();
      t_Clusters :: iterator i_old_list = old_cluster_list.begin();
      for (; i_old_list != i_old_list_last ; ++i_old_list )
      {
        for (types::t_int op=0; op<point_op.get_size(); op++)
        {
          // initialize a new cluster to object pointed by i_old_list
          Cluster transfo_cluster( *i_old_list );
          
          // transforms cluster according to symmetry group
          transfo_cluster.apply_symmetry(point_op(op),trans(op));
  
          // checks wether transformed cluster is in Lamarck::clusters
          t_Clusters :: iterator i_cluster  = _out.begin();
          t_Clusters :: iterator i_last = _out.end();
          for ( ; i_cluster != i_last ; ++i_cluster)
            if ( transfo_cluster->equivalent_mod_cell(*i_cluster, inv_cell) ) 
              break;
  
          // if it isn't, adds cluster to clusters
          if (i_cluster == i_last ) 
            _out.push_back( transfo_cluster );
        }
      }
    }
  }
  
} // namespace CE

