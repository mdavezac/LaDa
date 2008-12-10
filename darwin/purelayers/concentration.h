//
//  Version: $Id$
//
#ifndef _TWOSITES_H_
#define _TWOSITES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>

#include "concentration.h"

namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {
      //! \brief All concentration related behaviors
      //! \details Works for a single, set cell-shape. The concentration of the two
      //!          sites can be fixed. They can also be linked through a
      //!          load-minimizing function (see X_vs_y).  Finally, they can be
      //!          independant. In cases where the concentration is fully or
      //!          partially constrained, this class is able to set the right
      //!          concentration of an Crystal::Structure instance, and of a
      //!          TwoSites::Object instance.  It works correctly with structure
      //!          for which the occupation of some sites are frozen.
      //! \xmlinput see X_vs_y.
      class Concentration : public X_vs_y
      {
        public:
          //! Concentration of the first site. Output variable for Concentration::get()
          types::t_real x;
          //! Concentration of the second site. Output variable for Concentration::get()
          types::t_real y;
          //! The number of sites in this cell-shape.
          types::t_unsigned N;
          //! The number of frozen first lattice sites
          types::t_int Nfreeze_x;
          //! The number of frozen second lattice sites
          types::t_int Nfreeze_y;
          //! \brief Record wether the components of a Object::container is a first or
          //!        second lattice site.
          std::vector<bool> sites;
 
        public:
          //! Constructor
          Concentration  () 
                        : X_vs_y(), x(0), y(0), N(0),
                          Nfreeze_x(0), Nfreeze_y(0) {}
          //! Copy Constructor
          Concentration   ( const Concentration &_conc)
                        : X_vs_y(_conc), x(_conc.x), y(_conc.y),
                          N(_conc.N), Nfreeze_x(_conc.Nfreeze_x), Nfreeze_y(_conc.Nfreeze_y),
                          sites(_conc.sites) {}
          //! Destructor
          ~Concentration() {}
 
          //! \brief Loads the required behavior from XML (eg constrained,
          //!        partially constrained...)
          bool Load( const TiXmlElement &_node );
 
          //! \brief Normalizes the site occupations as given by the \b k-vectors. 
          //! \details A "complex" real-space occupation is computed from the
          //!          k-space intensities. For first-lattice sites, normalized
          //!          site-occupations are set to +/- 1 depending on wether the
          //!          complex real-space occupation is on left or right complex
          //!          half-plane. For second-lattice sites, normalized
          //!          site-occupations are set to +/- 1 depending on wether the
          //!          complex real-space occupation is on upper or lower complex
          //!          half-plane. These behaviors are in-line with the encoding
          //!          expected by the fourier transform.  If the concentration is
          //!          fixed, those sites for which the real value of the complex
          //!          occupation  are closest to zero are flipped first.
          //! \see GA::Krossover, GA::KRandom, GA::KMutation.
          void operator()( Crystal::Structure &_str );
          //! \brief Sets the concentration of \a _obj  if the concentration
          //!        is fully or partially constrained
          template< class T_OBJECT > void operator()( T_OBJECT &_obj );
          //! \brief Computes the concentration of \a _str and stores the result in
          //!        Concentration::x and Concentration::y
          void get( const Crystal::Structure &_str);
          //! \brief Computes the concentration of \a _obj and stores the result in
          //!        Concentration::x and Concentration::y
          template< class T_OBJECT > void get( const T_OBJECT &_obj );
          //! \brief Computes the number of sites in \a _str for which the
          //!        occupation is fixed
          void setfrozen ( const Crystal::Structure &_str );
 
        protected:
          //! \brief Actually does the job of normalizing the site occupations for
          //!        Concentration::operator()(Crystal::Structure&)
          void normalize( Crystal::Structure &_str, const types::t_int _site, 
                          types::t_real _tochange);
 
      };

    } // namespace PureLayers

    //! Policies for keeping track of functional evaluations in objects.
    namespace Keepers
    {
      struct ConcTwo
      {
        friend class boost::serialization::access;
        //! The concentration of the first site.
        types::t_real x;
        //! The concentration of the second site.
        types::t_real y;
        //! Constructor.
        ConcTwo() : x(-2), y(-2) {}
        //! Copy Constructor.
        ConcTwo( const ConcTwo& _c ) : x(_c.x), y(_c.y) {}
        //! Loads from attributes of \a _node.
        bool Load( const TiXmlElement &_node );
        //! Saves as attributes of \a _node.
        bool Save( TiXmlElement &_node ) const;
        private:
          //! Serializes concentration class.
          template<class Archive>
            void serialize(Archive & _ar, const unsigned int _version)
              { _ar & x; _ar & y; }
      };
    } // namespace Keepers
  } // namespace GA
} // namespace LaDa


#endif // _TWOSITES_OBJECT_H_
