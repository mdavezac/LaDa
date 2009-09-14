//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/lattice.h>

#include "find_all_cells.h"
#include "remove_smith_equivalents.h"
#include "exceptions.h"
#include "translation.h"
#include "label_exchange.h"


namespace LaDa
{
  namespace enumeration
  {
    //! Turns off translational equivalents in a SmithGroup.
    void remove_smith_equivalents( Database &_database,
                                   SmithGroup const &_sg, 
                                   Crystal::Lattice const &_lattice )
    {
      size_t const nsites( _lattice.sites.size() );
      size_t const card( _sg.smith(0)*_sg.smith(1)*_sg.smith(2) * nsites );
      size_t const nflavors( count_flavors(_lattice) );
      check_size(card, nflavors);
      if( nsites < 2 ) return;
      if( _database.size() != std::pow( count_flavors(_lattice), card ) )
        BOOST_THROW_EXCEPTION( incorrect_database() );
      
      boost::shared_ptr< std::vector<Translation> > translations = create_translations( _sg.smith, nsites );
      boost::shared_ptr< FlavorBase > flavorbase = create_flavor_base( card, nsites );
 
      typedef std::vector<Translation> :: const_iterator t_cit; 
      LabelExchange label_exchange(card, nflavors);
      t_cit i_begin = translations->begin();
      t_cit i_end = translations->end();
      for( t_uint x(0), max(std::pow<t_uint>(nflavors, card)); x < max; ++x )
      {
        if( not _database[x] ) continue;
        
        // removes label exchange
        while(++label_exchange) // bypass identity when entering loop.
        {
          t_uint const t( label_exchange(x, *flavorbase) );
#         ifdef LADA_DEBUG 
            if( t < x ) BOOST_THROW_EXCEPTION( internal() << error_string("t < x") ); 
#         endif
          _database[t] = false;
        }

        // loops over possible translations.
        for(t_cit i_trans(i_begin); i_trans != i_end; ++i_trans)
        {
          // removes translational and superperiodic equivalents, as well as their label exchange.
          t_uint t( (*i_trans)(x, *flavorbase) );
          // since the null translation is not included, t == x <=> superperiodic. 
#         ifdef LADA_DEBUG 
            if( t < x ) throw internal() << error_string("t < x"); 
#         endif
          _database[t] = false;
          
          // removes label exchange
          while(++label_exchange) // bypass identity when entering loop.
          {
            t_uint const t2( label_exchange(t, *flavorbase) );
#           ifdef LADA_DEBUG 
              if( t2 < x ) BOOST_THROW_EXCEPTION( internal() << error_string("t2 < x") ); 
#           endif
            _database[t2] = false;
          } // loop over label exchange.
        } // loop over pure translations.
      } // loop over x
    }
  }
}

