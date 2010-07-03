#ifndef LADA_CE_MLCLUSTER_H
#define LADA_CE_MLCLUSTER_H

#include "LaDaConfig.h"

#include <vector>
#include <iostream>

#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>

#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/symmetry_operator.h>
#include <crystal/smith.h>
#include <opt/propagate_std_vector.h>

namespace LaDa
{
  //! \brief Contains most everything %Cluster Expansion and structure related.
  //! \todo move structure related stuff to Physics
  //! \todo Rename Crystal to CE
  namespace CE
  {
    //! Defines a Multi-Lattice Cluster.
    class MLCluster 
    {
      friend class boost::serialization::access;
      
      public:
        //! Holds spin position of a cluster, relative to origin.
        struct Spin 
        {
          size_t site; //!< site index.
          math::rVector3d pos; //!< Position.
        };
        //! Type of the cluster spins.
        typedef std::vector<Spin> t_Spins;
        //! origin index.
        Spin origin;

        //! Constructor
        MLCluster () {}
        //! Copy Constructor
        MLCluster   (const MLCluster &_c)
                  : origin(_c.origin), spins_(_c.spins_) {}
        //! Destructor
        ~MLCluster () {}

        //! \brief Applies the point symmetry operation \a _op and the translation 
        //!        \a _trans to the cluster.
        void apply_symmetry(Crystal::SymmetryOperator const &_op);

        //! Returns the number of spins in the cluster (eg order of the interaction)
        types::t_int order() const { return size()+1; };

        //! Load a cluster from XML.
        bool load(const TiXmlElement &_node);
        //! True if involves same positions.
        bool operator==( MLCluster const & _c ) const;
        //! True if involves same positions.
        bool operator!=( MLCluster const & _c ) const
          { return not operator==(_c); }

        //! Computes pis for a given structure.
        types::t_real operator()( Crystal::Structure const &_str, 
                                  std::vector< std::vector<size_t> > const &_map,
                                  Crystal::t_SmithTransform const &_transform ) const;
        
        //! Computes pis for a given atom of a given structure.
        types::t_real operator()( Crystal::Structure::t_Atom const &_atom,
                                  Crystal::Structure const &_str, 
                                  std::vector< std::vector<size_t> > const &_map,
                                  Crystal::t_SmithTransform const &_transform ) const;
        // propagates (some)_std::vector members to this class.
        LADA_PROPAGATE_STDVECTOR(t_Spins, spins_);

      private:
        //! Vectors linking sites which interact through this cluster.
        t_Spins spins_;

        //! Serializes a cluster.
        template<class Archive> void serialize( Archive & _ar, const unsigned int _version)
          { _ar & origin; _ar & spins_; }
    };

    //! Compares two spins.
    bool operator==(MLCluster::Spin const &_a, MLCluster::Spin const &_b);
    //! Compares two spins.
    inline bool operator!=(MLCluster::Spin const &_a, MLCluster::Spin const &_b)
      { return not( _a == _b ); }

    //! Dumps a spin to as stream.
    std::ostream &operator<<( std::ostream &_sstr, const MLCluster::Spin &_spin );
    //! Dumps a cluster to as stream.
    std::ostream &operator<<( std::ostream &_sstr, const MLCluster &_cls );

  } // namespace CE

} // namespace LaDa

namespace boost {
  namespace serialization {

    //! Serializes a cluster spin.
    template<class Archive>
    void serialize(Archive & ar, LaDa::CE::MLCluster::Spin & g, const unsigned int version)
     { ar & g.site; ar & g.pos; }
  }
}
#endif
