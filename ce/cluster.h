//
//  Version: $Id$
//
#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <boost/filesystem/path.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
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
    // Forward declarations.
    //! \cond
      template<class T_HARMONIC> class Builder;  
      class Cluster;  
    //! \endcond


    //! Reads clusters from a NREL Jfile.
    void read_clusters( const Crystal::Lattice &_lat, 
                        const boost::filesystem::path &_path, 
                        std::vector< std::vector< Cluster > > &_out,
                        const std::string & _genes = "" );

    //! Reads a single cluster.
    bool read_cluster( const Crystal::Lattice &_lat, 
                       std::istream & _sstr,
                       Cluster &_out );
    
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
    class Cluster 
    {
      friend class boost::serialization::access;
      
      public:
        //! Vectors linking sites which interact through this cluster.
        std::vector<Eigen::Vector3d> vectors;

        //! Interaction energy of the cluster.
        types::t_real eci;

        //! Constructor
        Cluster () { eci = 0; };
        //! Copy Constructor
        Cluster (const Cluster &_c) : vectors( _c.vectors ), eci( _c.eci) {}
        //! Destructor
        ~Cluster () {}

        //! Sets the cluster interaction energy.
        types::t_real set_eci(const types::t_real _eci) { eci = _eci; return eci; };
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
        bool are_periodic_images( Cluster &_cluster, const Eigen::Matrix3d &_icell);

        //! Returns a constant reference to  vector \a _i of the cluster.
        const Eigen::Vector3d& access(types::t_int _i) const { return vectors[_i]; }
        //! Returns a reference to  vector \a _i of the cluster.
        Eigen::Vector3d& access(types::t_int _i) { return vectors[_i]; }
        //! Returns a constant reference to  vector \a _i of the cluster.
        const Eigen::Vector3d& operator[](types::t_int _i) const { return vectors[_i]; }
        //! Returns a reference to  vector \a _i of the cluster.
        Eigen::Vector3d& operator[](types::t_int _i) { return vectors[_i]; }

        //! Returns the number of vectors in the cluster (eg order of the interaction)
        types::t_int size() const { return vectors.size(); };

        //! Dumps the cluster to a string.
        void print_out( std::ostream &stream ) const; 
        //! Load a cluster from XML.
        bool Load(const TiXmlElement &_node);
        //! Returns  a constant reference to the vectors.
        const std::vector< Eigen::Vector3d >& Vectors() const { return vectors; }
        //! Returns  a reference to the vectors.
        std::vector< Eigen::Vector3d >& Vectors() { return vectors; }
        //! True if involves same positions.
        bool operator==( Cluster const & _c ) const;
        //! True if involves same positions.
        bool operator!=( Cluster const & _c ) const
          { return not operator==(_c); }
      private:
        //! Serializes a cluster.
        template<class Archive> void serialize( Archive & _ar,
                                                const unsigned int _version)
          { _ar & eci; _ar & vectors; } 
    };

    //! A single class of equivalent clusters.
    typedef std::vector<Cluster> t_Clusters;
    //! A vector of classes of equivalent clusters.
    typedef std::vector<t_Clusters> t_ClusterClasses;

    inline std::ostream &operator<<( std::ostream &_sstr, const Cluster &_cl )
      { _cl.print_out( _sstr ); return _sstr; }

    //! Adds clusters equivalent by symmetry.
    void  add_equivalent_clusters( const Crystal::Lattice &_lat, 
                                   std::vector< Cluster > &_out  );
    

  } // namespace CE

} // namespace LaDa

#endif
