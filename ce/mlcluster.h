//
//  Version: $Id$
//
#ifndef LADA_CE_MLCLUSTER_H
#define LADA_CE_MLCLUSTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <boost/filesystem/path.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <atat/machdep.h>
#include <atat/vectmac.h>
#include <mpi/mpi_object.h>

#include <crystal/lattice.h>
#include <crystal/structure.h>

namespace LaDa
{
  //! \brief Contains most everything %Cluster Expansion and structure related.
  //! \todo move structure related stuff to Physics
  //! \todo Rename Crystal to CE
  namespace CE
  {
    /** \brief Defines a cluster of a %Cluster Expansion %Functional
     *  \details A Cluster consists of an interaction energy Cluster::eci
     *           (\f$J\f$) and of \e N vectors Cluster::vectors
     *           (\f$\overrightarrow{v}_i\f$) linking a  site
     *           \f$S_{\overrightarrow{0}}\f$ with other sites. The resulting
     *           energy \f$E(S_{\overrightarrow{0}})\f$ of site
     *           \f$S_{\overrightarrow{0}}\f$ is simply,
     *           \f[
     *                E(S_0) = J\sum_{i=0}^N S_{\overrightarrow{0}} \cdot
     *                                       S_{\overrightarrow{0}+\overrightarrow{v}_i}.
     *           \f].
     *           A figure (in NREL parlance) is comprised of all clusters which
     *           are equivalent through the point-symmetry operations of the lattice.
     *  \pre As for any %Cluster Expansion, a fixed lattice is implied throughout.
     *  \todo add a pointer to the lattice so that things are more self container?
     *  \todo Make this a collection of equivalent clusters.
     */ 
    class MLCluster 
    {
      friend class boost::serialization::access;
      
      public:
        //! Type of the vector of positions.
        typedef std::vector<atat::rVector3d> t_Positions;
        //! Type of the vector of site indices.
        typedef std::vector<site> t_Sites;
        //! origin index.
        size_t origin;
        //! Indices of the sites.
        t_Sites sites;
        //! Vectors linking sites which interact through this cluster.
        t_Positions vectors;
        //! is null cluster.
        bool is_null;

        //! Constructor
        Cluster () : origin(0), is_null(true) {};
        //! Copy Constructor
        Cluster   (const Cluster &_c)
                : origin(_c.origins), sites(_c.sites), 
                  vectors(_c.vectors ), eci(_c.eci), is_null(_c.is_null) {}
        //! Destructor
        ~Cluster () {}

        //! \brief Applies the point symmetry operation \a _op and the translation 
        //!        \a _trans to the cluster.
        void apply_symmetry(Crystal::SymmetryOperator const &_op);
        //! \brief Checks wether this instance and \a _cluster are equivalent in
        //!        cell-shape specialization defined by \a _icell.
        //! \details Cluster Expansion functionals are most usefull when
        //!          specialized to a cell-shape. When this happens, a cluster
        //!          can be equivalent to another \e via translational
        //!          symmetries. This %function returns true when this instance
        //!          and \a _cluster are equivalent in the cell-shape defined by
        //!          the inverse of \a _icell.
        //! \param _cluster some cluster
        //! \param _icell if \f$\Omega\f$ is the cell-shape, then \a _icell = \f$\Omega^{-1}\f$.

        //! Returns the number of vectors in the cluster (eg order of the interaction)
        types::t_int size() const { return vectors.size(); };

        //! Load a cluster from XML.
        bool Load(const TiXmlElement &_node);
        //! True if involves same positions.
        bool operator==( Cluster const & _c ) const;
        //! True if involves same positions.
        bool operator!=( Cluster const & _c ) const
          { return not operator==(_c); }
      private:
        //! Serializes a cluster.
        template<class Archive> void serialize( Archive & _ar,
                                                const unsigned int _version)
          { _ar & origin; _ar & sites; _ar & vectors; _ar & eci; } 
    };

    //! A single class of equivalent clusters.
    typedef std::vector<MLCluster> t_MLClusters;
    //! A vector of classes of equivalent clusters.
    typedef std::vector<t_MLClusters> t_MLClusterClasses;

    std::ostream &operator<<( std::ostream &_sstr, const Cluster &_cl );
  } // namespace CE

} // namespace LaDa

#endif
