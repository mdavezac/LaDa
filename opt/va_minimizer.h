//
//  Version: $Id$
//
#ifndef _MINIMIZER_VA_H_
#define _MINIMIZER_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/utils/eoRNG.h>

#include <opt/types.h>

#include "opt_minimize_base.h"

namespace minimizer
{
  template< class T_FUNCTIONAL > std::string print( const T_FUNCTIONAL &_func )
  {
    typename T_FUNCTIONAL :: const_iterator  i_var = _func.begin();
    typename T_FUNCTIONAL :: const_iterator  i_var_end = _func.end();
    std::string str = " Vars: ";
    for(; i_var != i_var_end; ++i_var )
      str.push_back( ( *i_var > 0 ? '1': '0' ) ); 
    return str;
  }

  template< class T_FUNCTIONAL >
  class SaveNoState
  { 
    public:
      SaveNoState() {}
      SaveNoState( const SaveNoState<T_FUNCTIONAL> &) {}
    public:
      void save() {};
      void reset() {};
  };

  template< class T_FUNCTIONAL, class T_SAVESTATE = SaveNoState<T_FUNCTIONAL> > 
  class VA : public Base<T_FUNCTIONAL>
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
      typedef T_SAVESTATE  t_SaveState;
      typedef typename t_Functional :: t_Type t_Type;
      typedef typename t_Functional :: t_Container :: iterator iterator;
      typedef typename t_Functional :: t_Container :: const_iterator const_iterator;

    protected:
      using Base<t_Functional> :: current_func;

    protected:
      std::vector<types::t_unsigned> directions;
      types::t_unsigned itermax;
      bool do_check_gradient;
      t_SaveState *save_state;

    public:
      VA   ( const TiXmlElement &_node )
         : itermax(0), do_check_gradient(true),save_state(NULL) { Load(_node); }
      VA   ( const TiXmlElement &_node, t_SaveState &_s )
         : itermax(0), do_check_gradient(true), save_state(&_s) { Load(_node); }
      VA   ( const VA<t_Functional, t_SaveState> &_c ) 
         : itermax(_c.itermax),
           do_check_gradient( _c.do_check_gradient),
           save_state(_c.s) {}

      virtual ~VA() {}


      virtual bool operator()()
      {
        // then gets problem size of object
        types::t_unsigned size = current_func->size();

        // creates a vector of ordered directions
        directions.clear();
        { // implements sgl iota
          // -- i.e. fills which_direction with increasing intergers
          for( types::t_unsigned i =0; i < size; ++i)
            directions.push_back(i);
        }

        // now for the simmulated annealing part
        types::t_unsigned count = 0; // number of calls
        std::vector<types::t_unsigned> :: iterator i_begin = directions.begin();
        std::vector<types::t_unsigned> :: iterator i_last = directions.end();
        std::vector<types::t_unsigned> :: iterator i_dir;
        bool is_not_converged;
        iterator i_var_begin = current_func->begin();
        iterator i_var;
        t_Type current_e, next_e;
  
        // evaluates first position
        current_e = current_func->evaluate();
        if( save_state ) save_state->save();

        i_var = current_func->end();
        i_var_begin = current_func->begin();


        do 
        {
          ++count;

          // shuffle directions
          types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& ) = &eo::random<types::t_unsigned>;
          std::random_shuffle(i_begin, i_last, ptr_to_rng );
          is_not_converged = false;

          for( i_dir = i_begin; i_dir != i_last; ++i_dir )
          {
            i_var = i_var_begin + *i_dir; // sets direction

            // we only want investigate if direction is not taboo
            *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
            bool do_investigate = not current_func->is_taboo();
            *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
              
            if ( do_investigate  ) 
            {
              types::t_real grad = current_func->evaluate_one_gradient( *i_dir );
              if (     do_check_gradient 
                   and  grad < 0 )
              {
                *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
                if( save_state ) save_state->save();
                next_e = current_func->evaluate();
                if ( current_e > next_e )
                { // yeah! gradient was right
                  is_not_converged = true;
                  current_e = next_e;
                }
                else // flips back -- gradient was wrong
                {
                  *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
                  if( save_state ) save_state->reset();
                }
              } // end of gradient sign check
            } // end of check for taboo
            
            // breaks out after itermax calls to functional
            // unless itermax is null
            if ( itermax and count >= itermax )  
              { is_not_converged = false; break; }

          } // loop over shuffled directions

        } while ( is_not_converged );

        return true;
      } // end of minimize()

      virtual bool Load(const TiXmlElement &_node )
      {
        const TiXmlElement *parent;
        std::string str;
      
        // This whole section tries to find a <Functional type="vff"> tag
        // in _node or its child
        str = _node.Value();
        if ( str.compare("Minimizer" ) ) // returns zero if strings match exactly
          parent = _node.FirstChildElement("Minimizer");
        else
          parent = &_node;
      
        
        while (parent)
        {
          if ( parent->Attribute( "type" )  )
          {
            str = parent->Attribute("type");
            if ( str.compare("VA") == 0)
            {
              do_check_gradient = true;
              break;
            }
            else if ( str.compare("SA") == 0 )
            {
              do_check_gradient = false;
              break;
            }
          }
          parent = parent->NextSiblingElement("Minimizer");
        }
        if ( not parent )
        {
          std::cerr << "Could not find an <Minimizer type=\"VA or SA\"> tag in input file" 
                    << std::endl;
          return false;
        } 

        if ( parent->Attribute("itermax") )
        {
          types::t_int n = 0;
          parent->Attribute( "itermax", &n );
          itermax = ( n > 0 ) ? (types::t_unsigned) std::abs(n) : 0;
        }

        return true;
      }

  };

  template <class T_FUNCTIONAL>
  class Beratan : public Base<T_FUNCTIONAL>
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
      typedef typename t_Functional :: t_Type t_Type;
      typedef typename t_Functional :: t_Container :: iterator iterator;
      typedef typename t_Functional :: t_Container :: const_iterator const_iterator;

    protected:
      using Base<t_Functional> :: current_func;

    protected:
      std::vector<types::t_unsigned> directions;
      types::t_unsigned itermax;
      bool do_check_gradient;

    public:
      Beratan   ( const TiXmlElement &_node )
              : itermax(0), do_check_gradient(true) { Load(_node); }
      Beratan   ( const VA<t_Functional> &_c ) 
              : itermax(_c.itermax),
                do_check_gradient( _c.do_check_gradient) {}

      virtual ~Beratan() {}

      virtual bool operator() ()
      {
        iterator i_var_begin = current_func->begin();
        types::t_unsigned N = current_func->size(); // number of calls
        types::t_unsigned count = 0; // number of calls
        bool is_not_converged;

        // now for the simmulated annealing part
        t_Type *gradient = new t_Type[N];
        if ( not gradient )
        {
          std::cerr << "Memory Allocation error in minimizer::Beratan"    
                    << std::endl;
          return false;
        }
        t_Type *i_grad;
        t_Type *i_grad_last = gradient + N;
  
        // evaluates first position
        t_Type current_e = current_func->evaluate_with_gradient( gradient );

        do 
        {
          ++count;
          is_not_converged = false;
          
          // computes all gradient directions

          // and finds non-taboo direction with minimum gradient
          types::t_unsigned n = 0;
          do
          {
             i_grad = gradient;
             t_Type *min_dir = std::min_element( i_grad, i_grad_last );
             // if all remaining gradients are positive
             // convergence has been achieved
             if ( *min_dir >= 0 ) { delete[] gradient; return true; }
             
             // checks for taboo
             iterator i_var = i_var_begin + (min_dir - gradient); // sets direction
             *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
             if( not current_func->is_taboo() )
             {
               t_Type next_e = current_func->evaluate();
               if( current_e > next_e )
               {
                 current_e = next_e;
                 is_not_converged = true; 
                 break; 
               }
             }
             *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
             *min_dir = 1.0;
             ++n;

          } while (n <= N);


          // breaks out after itermax calls to functional
          // unless itermax is null
          if ( itermax and count >= itermax )  
            { is_not_converged = false; break; }

          current_func->evaluate_gradient( gradient );

        } while ( is_not_converged );


        delete[] gradient;
        return true;
      } // end of minimize()

      virtual bool Load(const TiXmlElement &_node )
      {
        const TiXmlElement *parent;
        std::string str;
      
        // This whole section tries to find a <Functional type="vff"> tag
        // in _node or its child
        str = _node.Value();
        if ( str.compare("Minimizer" ) != 0 )
          parent = _node.FirstChildElement("Minimizer");
        else
          parent = &_node;
      
        
        while (parent)
        {
          if ( parent->Attribute( "type" )  )
          {
            str = parent->Attribute("type");
            if ( str.compare("Beratan") == 0 )
              break;
          }
          parent = parent->NextSiblingElement("Minimizer");
        }
        if ( not parent )
        {
          std::cerr << "Could not find an <Minimizer type=\"VA or SA\"> tag in input file" 
                    << std::endl;
          return false;
        } 

        if ( parent->Attribute("itermax") )
        {
          types::t_int n = 0;
          parent->Attribute( "itermax", &n );
          itermax = ( n > 0 ) ? (types::t_unsigned) std::abs(n) : 0;
        }

        return true;
      }

  };
}
#endif
