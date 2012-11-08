#include <iostream>
#include <sstream>
#include <string>

#include <misc/types.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>
#include <math/misc.h>

#include "functional_builder.h"

namespace LaDa
{

  namespace CE 
  {
     template< class T_HARMONIC > Crystal::Lattice*
       Builder<T_HARMONIC>::lattice     = NULL;
//    template< class T_HARMONIC > typename Builder<T_HARMONIC> :: t_Clusters*
//      Builder<T_HARMONIC>::clusters    = NULL;
     template< class T_HARMONIC > typename Builder<T_HARMONIC> :: t_CS*
       Builder<T_HARMONIC>::harmonics   = NULL;

     template< class T_HARMONIC >
     Builder<T_HARMONIC> :: ~Builder<T_HARMONIC>()
     {
       if ( lattice )   delete lattice; 
       if ( clusters )  delete clusters; 
       if ( harmonics ) delete harmonics; 
       lattice = NULL; clusters = NULL; harmonics = NULL;
     }


     template< class T_HARMONIC >
     bool Builder<T_HARMONIC> :: Load (const TiXmlElement &_node)
     {
       const TiXmlElement *child, *parent;
       math::rVector3d vec;
       
       parent = opt::find_node( _node, "Functional", "type", "CE" );
       LADA_DO_NASSERT( not parent,
                   "Could not find an <Functional type=\"CE\"> tag in input file.\n" )

       // creates static quantities
       LADA_TRY_CODE(
         lattice    = new Crystal::Lattice;
         clusters   = new t_Clusters;
         harmonics  = new t_CS;,
         " Memory allocation failure in Builder :: Load (...) \n" 
       )

       // then loads lattice
       child = parent->FirstChildElement( "Lattice" );
       if ( not child )
       {
         const TiXmlNode *grandparent = parent->Parent();
         if ( parent )
           child = grandparent->FirstChildElement( "Lattice" );
         LADA_DO_NASSERT( not child, "Could not find Lattice in input.\n" )
       }
       LADA_DO_NASSERT( not lattice->Load(*child),
                   "Error while reading Lattice from input.\n" )
       lattice->find_space_group();
       Crystal :: Structure :: lattice = lattice;
       
       // then load clusters
       child = parent->FirstChildElement( "Clusters" );
       LADA_DO_NASSERT( not child,
                   "Could not find Clusters in input.\n" )
       child = child->FirstChildElement( "Cluster" );
       t_Cluster cluster;
       for (  ; child; child = child->NextSiblingElement("Cluster") )
       {
         LADA_DO_NASSERT( not cluster.Load( *child ),
                     "Error while loading cluster from input.\n" )
         clusters->push_back(cluster);
       }
       
       // reads in Constituent Strain (coefficients are static)
       child = parent->FirstChildElement( "CS" );
       LADA_DO_NASSERT( not child,
                   "Could not find CS in input.\n" )
       LADA_DO_NASSERT( not harmonics->Load_Harmonics( *child ),
                   "Error while loading harmonics from input.\n" )

       return true;
     }


     // converts the clusters into a polynomial forme
     // converts each cluster for each atom to a monome
     // One function per structure
     // each function is made up of a polynomial and constituent strain
     // part
     template< class T_HARMONIC >
     bool Builder<T_HARMONIC> :: generate_functional( const Crystal::Structure &str,
                                                     t_VA_Functional * const functional ) const
     {
       std::pair<t_Chemical*, t_CS*> pair( generate_functional( str ) );
       if( (not pair.first) or (not pair.second) ) 
       {
         if( pair.first ) delete pair.first;
         if( pair.second ) delete pair.second;
         return false;
       }
       functional->set_functional2( pair.second );
       functional->set_functional1( pair.first );
   
       return true;
     }

     template< class T_HARMONIC >
     std::pair< typename Builder<T_HARMONIC> :: t_Chemical*,
                typename Builder<T_HARMONIC> :: t_CS* >
       Builder<T_HARMONIC> :: generate_functional( const Crystal::Structure &str ) const
       {
         t_Chemical *polynome;
    
         // finally, creates polynomials
    
         math::rMatrix3d inv_cell = str.cell.inverse();
         polynome = new t_Chemical();
         Crystal :: Structure :: t_Atoms :: const_iterator i_atom = str.atoms.begin();
         Crystal :: Structure :: t_Atoms :: const_iterator i_atom_last = str.atoms.end();
    
         for(; i_atom != i_atom_last; ++i_atom) // loop over atoms
         {
           t_Clusters :: iterator i_cluster = clusters->begin();
           t_Clusters :: iterator i_cluster_last = clusters->end();
           math::rVector3d atom_pos = i_atom->pos;
           
           for( ; i_cluster != i_cluster_last; ++i_cluster ) // loop over clusters
           {
             typedef std::vector<math::rVector3d> :: iterator vec_iterator;
             vec_iterator i_cpos_begin = i_cluster->vectors.begin();
             vec_iterator i_cpos_center = i_cluster->vectors.begin();
             vec_iterator i_cpos_last = i_cluster->vectors.end();
             vec_iterator i_cpos;
             if ( i_cluster->vectors.size() == 0 )
             {
               i_cpos_center = i_cpos_last;
               polynome->add( function::Monome<>(i_cluster->eci) ); 
             }
    
             for (; i_cpos_center != i_cpos_last; ++i_cpos_center ) 
             {   
               // sets up a monome with the correct coefficient
               function::Monome<> monome(   i_cluster->eci
                                          / ( (types::t_real) i_cluster->vectors.size() ) );
    
               i_cpos = i_cpos_begin;
               for (; i_cpos != i_cpos_last; ++i_cpos ) // loop over cluster points
               {
                 math::rVector3d shift = atom_pos - *i_cpos_center;
                 
                 if ( not math::is_integer( (lattice->cell.inverse()*shift).eval() ) ) continue;
                 
                 // finds atom to which "point" is equivalent
                 Crystal::Structure::t_Atoms::const_iterator i_equiv = str.atoms.begin();
                 size_t index(0);
                 for (; i_equiv != i_atom_last; ++i_equiv, ++index)  
                   if ( math::are_periodic_images(*i_cpos + shift, i_equiv->pos, inv_cell) ) 
                     break;
                 
                 LADA_NASSERT( i_equiv == i_atom_last, 
                             "Could not find equivalent of atom "
                           << (*i_cpos + shift) << " in Lamarck::generate_functionals\n" )
                 monome.add_term(index, true);
                 
               }  // end of loop over cluster points
    
               polynome->add(monome);
               
             } // end of rotation
    
             
           }  // end of loop over clusters
         }  // end of loop over atoms 
    
         (*polynome) *= 1.0 /( (types::t_real) str.atoms.size() );
    
         // now computes constituent strain 
         return std::pair< t_Chemical*, t_CS* >( polynome, new t_CS(str) );
       }
  } // namespace CE
} // namespace LaDa


