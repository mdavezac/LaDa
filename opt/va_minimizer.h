//
//  Version: $Id$
//
#ifndef _MINIMIZER_VA_H_
#define _MINIMIZER_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/utils/eoRNG.h>

#include "types.h"
#include "fuzzy.h"

#include "minimize_base.h"

namespace Minimizer
{
  //! \brief Dummy save-state Functor. 
  //! \details Save-state functors should be capable of, well, saving the
  //!          current state within the minimization process so that the
  //!          minimizer can go back to it later on. This implements a dummy
  //!          stack for when no save-state functor  is provided. It can be
  //!          used by minimizers for which saving a state is optional. 
  //! \note This functor provides all necessary behaviors for a save-state
  //!       functor. Copy it!
  template< class T_FUNCTIONAL >
  class SaveNoState
  { 
    public:
      //! Constructor
      SaveNoState() {}
      //! Copy Constructor
      SaveNoState( const SaveNoState<T_FUNCTIONAL> &) {}
    public:
      //! Dummy save.
      void save() {};
      //! Dummy reset.
      void reset() {};
  };

  //! \brief Virtual Atom minimizer as described <A
  //!        HREF="http://dx.doi.org/10.1088/0953-8984/19/40/402201"> here  </A>.
  //! \details This is a minimizer for discrete functionals. It checks the
  //!          gradient in random directions. If the gradient is negative, then
  //!          it computes the value of the functional in one step in that
  //!          direction. If the resulting value is indeed lower, it proceeds
  //!          from there. Otherwise, it goes to check another direction.
  //!          Convergence is deemed achive when all directions have been
  //!          checked without finding a better neighbor.
  //!           
  //!          The directions to compute the gradient for are chosen for a
  //!          complete of list of randomly shuffled directions.
  //!          
  //!          In addition, before considering a direction, the minimizer will
  //!          check that the neighboring point in that direction is not taboo.
  //! 
  //!          On input the options are:
  //!             - Maximum number of iterations before checking out.
  //!             - Wether or not to check the gradient in a direction before
  //!               computing the value of that neighbor.
  //!             - Wether to use a save-state functor. This "option" is set by
  //!               including a save-state functor directly in the arguments of
  //!               the constructor.
  //!             .
  //! \xmlinput 
  //!   \code 
  //    <Minimizer type="VA" itermax=100/>
  //!   \endcode.
  //!   The %VA minizer does check the gradients before moving to a neighboring
  //!   point. The %SA minizer does not. The attribute \a itermax specifies the
  //!   maximum number of iterations to do before checking out. 
  template< class T_FUNCTIONAL, class T_SAVESTATE = SaveNoState<T_FUNCTIONAL> > 
  class VA : public Base<T_FUNCTIONAL>
  {
    protected:
      //! The type of the base class.
      typedef Base<T_FUNCTIONAL>  t_Base;
    public:
      //! The type of the functional.
      typedef T_FUNCTIONAL t_Functional;
      //! The save-state functor.
      typedef T_SAVESTATE  t_SaveState;
      //! See function::Base::t_Type
      typedef typename t_Functional :: t_Type t_Type;

    protected:
      using Base<t_Functional> :: current_func;

    protected:
      //! The randomly shuffled list of directions.
      std::vector<types::t_unsigned> directions;
      //! The maximum number of iterations
      types::t_unsigned itermax;
      //! Wether or not to check the gradient
      bool do_check_gradient;
      //! A pointer to the save-state functor.
      t_SaveState *save_state;
      //! A pointer to another save-state functor.
      t_SaveState *save_new_state;

    public:
      //! Constructor and Initializer
      VA   ( const TiXmlElement &_node )
        //! \cond
         : itermax(0), do_check_gradient(true), 
           save_state(NULL), save_new_state(NULL) { Load(_node); }
        //! \endcond
        
      //! Constructor and Initializer
      VA   ( const TiXmlElement &_node, t_SaveState &_s )
        //! \cond
         : itermax(0), do_check_gradient(true),
           save_state( &_s ), save_new_state(NULL) { Load(_node); }
        //! \endcond

      //! Copy Constructor.
      VA   ( const VA<t_Functional, t_SaveState> &_c ) 
        //! \cond
         : itermax(_c.itermax),
           do_check_gradient( _c.do_check_gradient),
           save_state(_c.save_state), save_new_state( _c.save_new_state) {}
        //! \endcond

      //! Destructor
      virtual ~VA() {}


      //! The minimization functor.
      bool operator()();

      //! \brief Loads parameters from XML.
      virtual bool Load(const TiXmlElement &_node );

  };

  //! \brief Virtual Atom minimizer as described in <A
  //!        HREF="http://dx.doi.org/10.1021/jp0646168"> Keinan, Hue, Beratan,
  //!        and Yang, J. Phys. Chem A \b 111, 146 (2006)  </A>.
  //! \details This is a minimizer for discrete functionals. It checks the
  //!          gradient in all directions from a starting point. It moves in
  //!          the direction with the lowest gradient and check the value at
  //!          that neighboring point. If the value is indeed lower, it
  //!          iterates from there.  Otherwise it chooses the direction with
  //!          the next most negative gradient.  Convergence is deemed achive
  //!          when all directions have been checked without finding a better
  //!          neighbor.
  //!           
  //!          In addition, before considering a direction, the minimizer will
  //!          check that the neighboring point in that direction is not taboo.
  //! 
  //!          On input the options are:
  //!             - Maximum number of iterations before checking out.
  //!             .
  //! \xmlinput 
  //!   \code 
  //    <Minimizer type="Beratan" itermax=100/>
  //!   \endcode.
  //!   The attribute \a itermax specifies the maximum number of iterations
  //!   to do before checking out. 
  template <class T_FUNCTIONAL>
  class Beratan : public Base<T_FUNCTIONAL>
  {
    public:
      //! The type of the functional.
      typedef T_FUNCTIONAL t_Functional;
      //! See function::Base::t_Type
      typedef typename t_Functional :: t_Type t_Type;

    protected:
      using Base<t_Functional> :: current_func;

    protected:
      //! The list of directions
      std::vector<types::t_unsigned> directions;
      //! The maximum number of iterations
      types::t_unsigned itermax;

    public:
      //! Constructor and Initializer.
      Beratan   ( const TiXmlElement &_node )
              : itermax(0) { Load(_node); }
      //! Copy Constructor.
      Beratan   ( const VA<t_Functional> &_c ) 
              : itermax(_c.itermax) {}

      //! Destructor.
      virtual ~Beratan() {}

      //! Minimization functor.
      virtual bool operator() ();

      //! Loads the parameters from XML.
      virtual bool Load(const TiXmlElement &_node );

  };




  template< class T_FUNCTIONAL, class T_SAVESTATE>
    bool VA<T_FUNCTIONAL, T_SAVESTATE> :: operator()()
    {
      // initializes object related stuff
      if( not current_func )  return false;
      if( not current_func->init() )  return false;

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
      typename t_Functional :: t_Container :: iterator i_var_begin = current_func->begin();
      typename t_Functional :: t_Container :: iterator i_var;
      t_Type current_e, next_e;
      types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& );
      ptr_to_rng = &eo::random<types::t_unsigned>;
  
      // evaluates first position
      current_e = current_func->evaluate();
      if( save_state ) 
      {
        save_new_state = new T_SAVESTATE( *save_state );
        save_state->save();
      }
      i_var_begin = current_func->begin();


      do 
      {
        ++count;

        // shuffle directions
        std::random_shuffle(i_begin, i_last, ptr_to_rng );
        is_not_converged = false;

        for( i_dir = i_begin; i_dir != i_last; ++i_dir )
        {
          i_var = i_var_begin + *i_dir; // sets direction

          // we only want investigate if direction is not taboo
          *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
          bool do_investigate = not current_func->is_taboo();
          *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
            
          if ( not do_investigate  ) continue;
          
          if ( do_check_gradient )
          { 
            types::t_real grad =   current_func->evaluate_one_gradient( *i_dir ) 
                                 * ( *i_var > t_Type(0) ? -1e0: 1e0 );
            
            if ( Fuzzy::geq<t_Type>(grad, 0) ) continue;
          }

          if( save_state ) save_state->save();
          *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
          current_func->invalidate();
          next_e = current_func->evaluate();
          // Saves the state and recomputes current_e (eg for moving target
          // optimizations, such as convex hulls. Indeed in the case of a
          // convex hull, evaluating next_e may have added a breakpoint to the
          // convex hull and changed the target function.
          if( save_state )
          {
            save_new_state->save();
            *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
            save_state->reset();
            current_e = current_func->evaluate();
            save_state->save();
            *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
            save_new_state->reset();
          }

          if ( Fuzzy::geq(current_e, next_e) )
          {
            is_not_converged = true;
            current_e = next_e;
            if( save_state ) save_state->save();
            continue;
          }

          // flips back -- gradient was wrong
          *i_var = ( *i_var > t_Type(0) ) ? t_Type(-1) : t_Type(1); // flips spins
          if( save_state ) save_state->reset();

        } // loop over shuffled directions

        // breaks out after itermax loops over directions
        if ( itermax and count >= itermax )  is_not_converged = false;

      } while ( is_not_converged );

      if( save_new_state )
      {
        delete save_new_state; 
        save_new_state = NULL;
      }

      return true;
    } // end of minimize()

  template< class T_FUNCTIONAL, class T_SAVESTATE>
    bool VA<T_FUNCTIONAL, T_SAVESTATE> :: Load(const TiXmlElement &_node )
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







  template <class T_FUNCTIONAL>
   bool Beratan<T_FUNCTIONAL> :: operator() ()
    {
      // initializes object related stuff
      if( not current_func )  return false;
      if( not current_func->init() )  return false;

      typename t_Functional :: t_Container :: iterator i_var_begin = current_func->begin();
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
           typename t_Functional :: t_Container :: iterator i_var;
           i_var = i_var_begin + (min_dir - gradient); // sets direction
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

  template <class T_FUNCTIONAL>
    bool Beratan<T_FUNCTIONAL> :: Load(const TiXmlElement &_node )
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
}
#endif
