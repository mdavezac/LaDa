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
#include <opt/ndim_iterator.h>
#include <atat/machdep.h>
#include <atat/vectmac.h>
#include <mpi/mpi_object.h>

#include <crystal/lattice.h>

//! \cond
namespace CE
{ 
  template<class T_HARMONIC> class Builder; 
  class Cluster;  
};
//! \endcond

//! \brief Contains most everything %Cluster Expansion and structure related.
//! \todo move structure related stuff to Physics
//! \todo Rename Crystal to CE
namespace CE
{
  //! \brief Finds all clusters which are at most the \a _max_neigh and with at
  //!        most \a _maxN spins.
  //! \return A vector of vectors of cluster. The inner vectors contain
  //!         symmetrically equivalent clusters. The new cluster classes are
  //!         added to the output vector without checking for duplicates.
  //!
  void find_all_clusters( const Crystal :: Lattice &_lat,
                          types::t_unsigned _max_neigh,
                          types::t_unsigned _maxN,
                          std::vector< std::vector<Cluster> > &_out );
  //! Reads clusters from a NREL Jfile.
  void read_clusters( const Crystal::Lattice &_lat, 
                      const boost::filesystem::path &_path, 
                      std::vector< std::vector< Cluster > > &_out );

  
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
   */ 
  class Cluster 
  {
    template<class T_HARMONIC> friend class Builder;
    friend void find_all_clusters( const Crystal :: Lattice &,
                                   types::t_unsigned,
                                   types::t_unsigned,
                                   std::vector< std::vector<Cluster> > &);
    friend void read_clusters( const Crystal::Lattice &_lat,
                               const boost::filesystem::path &_path, 
                               std::vector< std::vector< Cluster > > &_out );
    
    protected:
      //! Vectors linking sites which interact through this cluster.
      std::vector<atat::rVector3d> vectors;

    public:
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
      void apply_symmetry(const atat::rMatrix3d &_op, const atat::rVector3d &_trans);
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
      bool equivalent_mod_cell( Cluster &_cluster, const atat::rMatrix3d &_icell);

      //! Returns a constant reference to  vector \a _i of the cluster.
      const atat::rVector3d& access(types::t_int _i) const { return vectors[_i]; }
      //! Returns a reference to  vector \a _i of the cluster.
      atat::rVector3d& access(types::t_int _i) { return vectors[_i]; }
      //! Returns a constant reference to  vector \a _i of the cluster.
      const atat::rVector3d& operator[](types::t_int _i) const { return vectors[_i]; }
      //! Returns a reference to  vector \a _i of the cluster.
      atat::rVector3d& operator[](types::t_int _i) { return vectors[_i]; }

      //! Returns the number of vectors in the cluster (eg order of the interaction)
      types::t_int size() const { return vectors.size(); };

      //! Dumps the cluster to a string.
      void print_out( std::ostream &stream ) const; 
      //! Load a cluster from XML.
      bool Load(const TiXmlElement &_node);
      //! Serializes a cluster.
      template<class Archive> void serialize( Archive & _ar,
                                              const unsigned int _version)
        { _ar & eci; _ar & vectors; } 
  };

  inline std::ostream &operator<<( std::ostream &_sstr, const Cluster &_cl )
    { _cl.print_out( _sstr ); return _sstr; }


 
} // namespace CE

#endif
