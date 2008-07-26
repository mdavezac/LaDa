//
//  Version: $Id$
//

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace CE
{
  template< class T_CLUSTERS, class T_PIS >
  void find_pis( const T_CLUSTERS &_clusters,
                 const Crystal::Structure & _str,
                 T_PIS &_pis )
  {
    typedef T_CLUSTERS  t_Clusters;
    typedef typename t_Clusters :: value_type t_EClusters;
    typedef T_PIS t_Pis;

    _pis.resize( _clusters.size() );

    atat::rMatrix3d inv_cell = !(~str.cell);
    Crystal :: Structure :: t_Atoms :: const_iterator i_atom = str.atoms.begin();
    Crystal :: Structure :: t_Atoms :: const_iterator i_atom_end = str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom) // loop over atoms
    {

      // loop over classes of clusters.
      t_EClusters :: const_iterator i_clusters = _clusters.begin();
      t_EClusters :: const_iterator i_clusters_end = _clusters.end();
      t_Pis :: iterator i_pi = pis.begin();
      for( ; i_clusters != i_clusters_end; ++i_clusters, ++i_pi ) 
      {
        *i_pi = 0;
        t_Clusters :: const_iterator i_cluster = i_clusters->begin();
        t_Clusters :: const_iterator i_cluster_end = i_clusters->end();
        // loop over equivalent clusters.
        for( ; i_cluster != i_cluster_end; ++i_cluster )
        {
          // null cluster
          if ( i_cluster->vectors.size() == 0 ) continue;

          // loop over cluster origin.
          typedef std::vector<atat::rVector3d> :: iterator vec_iterator;
          vec_iterator i_cpos_begin = i_cluster->vectors.begin();
          vec_iterator i_cpos_center = i_cluster->vectors.begin();
          vec_iterator i_cpos_end = i_cluster->vectors.end();
          vec_iterator i_cpos;
          for (; i_cpos_center != i_cpos_end; ++i_cpos_center ) 
          {   
            // loop over lattice positions in cluster.
            types::t_int result(1);
            for ( i_cpos = i_cpos_begin; i_cpos != i_cpos_end; ++i_cpos )
            {
              atat::rVector3d shift = i_atom->pos - *i_cpos_center;
              
              if (is_int( (!lattice->cell)*shift)) 
              {
                // finds atom to which lattice site is equivalent
                Crystal::Structure::t_Atoms::const_iterator
                  i_equiv = _str.atoms.begin()
                for (; i_equiv != i_atom_end; ++i_equiv)  
                  if ( atat::equivalent_mod_cell(*i_cpos + shift,
                                                 i_equiv->pos,inv_cell) ) 
                    break;
                
                __ASSERT( index == str.atoms.size(),
                            "Could not find equivalent of atom "
                          << (*i_cpos + shift)
                          << " in Lamarck::generate_functionals\n" )
                result *= i_equiv->type > 0 ? 1: -1;
              }

            }  // end of loop over cluster points
            *i_pi += result;
        
          } // end of rotation
        
        }  // end of loop over equivalent clusters
      } // end of loop over cluster classes.
    }  // end of loop over atoms 
  
    (*polynome) *= 1.0 /( (types::t_real) str.atoms.size() );
  
    // now computes constituent strain 
    return std::pair< t_Chemical*, t_CS* >( polynome, new t_CS(str) );
  }

  template< class T_CLUSTERS, class T_PIS >
  void find_pis( const T_CLUSTERS &_clusters,
                 const std::vector< Crystal::Structure > & _str;
                 T_PIS &_pis )
  {
    namespace bl = boost::lambda;
    _pis.resize( _str.size() );
    void (*ptr_func)( const T_CLUSTERS &_clusters,
                      const Crystal::Structure & _str;
                      T_PIS &_pis ) = &find_pis<T_CLUSTERS, T_PIS>;
    opt::concurrent_loop
    (
      _str.begin(), _str.end(), _pis.begin(),
      bl::bind( ptr_func, bl::constant( _clusters ), bl::_1, bl::_2 )
    );
  }

  template< template<class> T_CONTAINER >
    void read_ce_structures( const std::string &_path,
                             T_CONTAINER<Crystal::Structure> &_structures )
    {
      __TRYBEGIN
      namespace fs = boost::filesystem;

      // First finds directory of LDAs.dat.
      size_t t = pstring.find_last_of("/");
      fs::path filename(_path);
      fs::path path( t != std::string::npos ? _path.substr(0,t): "."  );

      // then starts reading file.
      std::ifstream ldas( fullpath.string().c_str(), std::ifstream::in );
      std::string line;
      while( std::getline( ldas, line ) )
      {
        const boost::regex re("^(\\s+)?(\\S+)\\s+(-?\\d+(\\.\\d+)?)");
        boost::match_results<std::string::const_iterator> what;
        if( not boost::regex_search( line, what, re ) ) continue;

        Crystal :: Structure structure;
        structure.name( what.str(2) );
        structure.energy( boost::lexical_cast<types::t_real>( what.str(3) ) );

        Crystal :: read_structure( structure, path / structure.name );
        _structures.push_back(structure)
      }

      __TRYEND(,"Error while parsing " << _path << " and structures.\n" )
    }



} // end of namespace CE

#endif 
} // end of namespace CE

