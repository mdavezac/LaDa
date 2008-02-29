//
//  Version: $Id$
//
#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <atat/machdep.h>
#include <atat/vectmac.h>
#include <mpi/mpi_object.h>

//! \cond
namespace VA_CE{ class Functional_Builder; };
namespace Ising_CE{ class Cluster; }
___DECLAREMPIOBJECT( Ising_CE::Cluster )
//! \endcond

//! \brief Contains most everything %Cluster Expansion and structure related.
//! \todo move structure related stuff to Physics
//! \todo Rename Ising_CE to CE
namespace Ising_CE
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
   */ 
  class Cluster 
  {
    //! \cond
    ___FRIENDMPIOBJECT(Ising_CE::Cluster)
    //! \endcond
    friend class VA_CE::Functional_Builder;
    
    protected:
      //! Vectors linking sites which interact through this cluster.
      std::vector<atat::rVector3d> vectors;
      //! Interaction energy of the cluster.
      types::t_real eci;

    public:
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
  };
  
} // namespace Ising_CE

#ifdef _MPI
#include <mpi/mpi_object.h>
#include <atat/serialize.h>
namespace mpi
{
#define ___OBJECTCODE \
   return     _this.serialize( _ob.eci ) \
          and _this.serialize_container( _ob.vectors );
#define ___TYPE__ Ising_CE::Cluster
  /** \ingroup MPI
  * \brief Serializes an Ising_CE::Cluster
  */
#include <mpi/serialize.impl.h>
#undef ___OBJECTCODE
}
#endif

#endif
