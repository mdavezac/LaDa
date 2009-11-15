//
//  Version: $Id$
//
#ifndef LADA_CE_MLCLUSTERS_H
#define LADA_CE_MLCLUSTERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/propagate_std_vector.h>
#include "mlcluster.h"

namespace LaDa
{
  namespace CE
  {
    //! Defines a vector of symmetry equivalent clusters.
    class MLClusters
    {
      friend class boost::serialization::access;
      
      public:
        //! Container of clusters.
        typedef std::vector<MLCluster> t_MLClusters;
        //! Interaction energy.
        types::t_real eci;

        //! Constructor
        MLClusters() : eci(0) {};
        //! Copy Constructor
        MLClusters   (const MLClusters &_c)
                   : clusters_(_c.clusters_), eci(_c.eci) {}
        //! Destructor
        ~MLClusters() {}

        //! Returns the number of spins in these clusters.
        types::t_int order() const { return empty() ? 0: back().order(); }

        //! Load a cluster from XML.
        bool load(const TiXmlElement &_node);
        //! True if this class contains _c.
        bool operator==(MLCluster const & _b) const;
        //! True if all clusters in _c are contained here, and vice-versa.
        bool operator==(MLClusters const & _b) const;
        //! False if this class contains _c.
        bool operator!=(MLCluster const & _b) const
          { return not operator==(_b); }
        //! False if all clusters in _c are contained here, and vice-versa.
        bool operator!=(MLClusters const & _b) const
          { return not operator==(_b); }

        //! Initializes a cluster class from a cluster.
        void init(MLCluster const &_cluster, types::t_real _eci=0e0);
        
        //! Computes pis for a given structure.
        types::t_real operator()( Crystal::Structure const &_str, 
                                  std::vector< std::vector<size_t> > const &_map,
                                  Crystal::t_SmithTransform const &_transform ) const;

        // propagates (some)_std::vector members to this class.
        LADA_PROPAGATE_STDVECTOR(t_MLClusters, clusters_);
        
      private:
        //! Container of clusters.
        t_MLClusters clusters_;
        //! Serializes a cluster.
        template<class Archive> void serialize( Archive & _ar, const unsigned int _version)
          { _ar & clusters_; _ar & eci; }
    };

    //! A vector of classes of equivalent clusters.
    typedef std::vector<MLClusters> t_MLClusterClasses;

    //! Dumps a cluster class to as stream.
    std::ostream &operator<<( std::ostream &_sstr, const MLClusters &_cls );

    //! Read clusters from NREL file
    boost::shared_ptr<t_MLClusterClasses> read_clusters( const Crystal::Lattice &_lat, 
                                                         const boost::filesystem::path &_path, 
                                                         const std::string & _genes = "" );
    //! Read clusters from NREL file
    boost::shared_ptr<t_MLClusterClasses> load_mlclusters( const boost::filesystem::path &_path, 
                                                           bool is_newformat = true );
    
  } // namespace CE

} // namespace LaDa
#endif
