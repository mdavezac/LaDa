#ifndef _OPT_MINIMIZE_LINEAR_H_
#define _OPT_MINIMIZE_LINEAR_H_

// linear minimizer: checks derivative at physical points, and jumps 
// Simulated Annealing: checks if other point is lower in E
// The two methods are equivalent if functional is linear with respect
// to each spin, although SA does require less calls.
// works only for physical points +/- 1 

#include <vector> 
#include <algorithm> 
#include <limits.h>

#include <opt/opt_minimize_base.h>
#include <analysis/analyze_code.h>
#include <opt/types.h>

namespace opt {

  template<typename t_Object> 
  class Minimize_Linear : public Minimize_Base< t_Object >
  {
    public: 
      static types::t_unsigned nb_guess;
      static types::t_unsigned good_guess;
      static types::t_unsigned bad_guess;
      static types::t_unsigned nb_evals;
      static types::t_unsigned nb_grad_evals;
      
    protected:
      using Minimize_Base<t_Object>::current_func;

    public:
      bool simulated_annealing; // T=0 simulated annealing
      types::t_unsigned max_calls;   // nb calls to evaluate
      
    public:
      inline Minimize_Linear( t_Object* _func ): Minimize_Base<t_Object>(_func),
                                               simulated_annealing(true), 
                                               max_calls(UINT_MAX){};
      inline Minimize_Linear(): Minimize_Base<t_Object>(),
                                simulated_annealing(true), 
                                max_calls(UINT_MAX) {};
      virtual ~Minimize_Linear() {};

      virtual bool minimize()
      {
        START_ANALYSIS( "Minimize :: minimize" )
        types::t_int pb_size = current_func->size();
        bool is_converged = false;
        typename t_Object::iterator i_var_begin = current_func->begin();
        typename t_Object::iterator i_var;
      
        std::vector<types::t_int> which_direction(pb_size);
        std::vector<types::t_int> :: iterator i_int, i_begin, i_last;
      
        { // implements sgl iota
          // -- i.e. fills which_direction with increasing intergers
          for( types::t_int i =0; i < pb_size; ++i)
            which_direction.push_back(i);
        }
      
        i_begin = which_direction.begin();
        i_last = which_direction.end();

        // convergence loop
        if ( simulated_annealing )
        {
          types::t_unsigned count = 0;
          do
          {
            // shuffles directions randomly;
            std::random_shuffle(i_begin, i_last);
        
            // sets convergence to true
            is_converged = true;
        
            // now iterates through directions and minimizes
            typename t_Object :: t_Type current_e = current_func->evaluate(); 
            typename t_Object :: t_Type other_e;
            ++count;
            ++nb_evals;

            for( i_int = i_begin; i_int != i_last; ++i_int )
            {
              i_var = i_var_begin + *i_int;
        
              // flips spin
              *i_var =  (*i_var > 0 ) ? -1.0 : 1.0;
              other_e = current_func->evaluate();
              ++count;
              ++nb_evals;
        
              // checks wether we should change value or not
              if ( current_e > other_e )
              {
                is_converged = false;
                current_e = other_e;
              }
              else
                *i_var =  (*i_var > 0 ) ? -1.0 : 1.0; // flips spin back
              if ( count >= max_calls )    // breaks out after n calls to functional
                { is_converged = true; break; }

            }
          }
          while ( !is_converged );
        } // gradient directed minimization
        else
        {
          types::t_unsigned count = 0;
          do
          {
            // shuffles directions randomly;
            std::random_shuffle(i_begin, i_last);
        
            // sets convergence to true
            is_converged = true;
        
            // now iterates through directions and minimizes
            typename t_Object :: t_Type current_e = current_func->evaluate(); 
            ++nb_evals;
            ++count;
            typename t_Object :: t_Type other_e;
            for( i_int = i_begin; i_int != i_last; ++i_int )
            {
              i_var = i_var_begin + *i_int;
        
              // computes derivative in direction i_int
              *i_var = ( ( *i_var > 0 ) ?  0.998: -0.998 ); 
              other_e = current_func->evaluate();
              ++count;
              ++nb_grad_evals;
        
              // checks wether we should change value or not
              if ( current_e > other_e )
              {
                *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // goes to "other" point
                other_e = current_func->evaluate();
                ++count;
                ++nb_evals;
                if ( current_e > other_e )
                {
                  is_converged = false;
                  current_e = other_e;
                  ++good_guess;
                }
                else
                {
                  *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // flips back
                  ++bad_guess;
                }
              }
              else 
                *i_var = ( *i_var > 0.0 ) ? 1.0 : -1.0; // goes to "other" point
            // else
            // {
            //   *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // goes to "other" point
            //   other_e = current_func->evaluate();
            //   if ( current_e <= other_e )
            //     ++good_guess;
            //   else 
            //     ++bad_guess;
            //   *i_var = ( *i_var > 0.0 ) ? -1.0 : 1.0; // goes to "other" point
            // }

              if ( count >= max_calls ) 
                { is_converged = true; break; }
            }
          }
          while ( !is_converged );
        } // end of derivative check

        END_ANALYSIS;
        return true;
      
      } // minimize
  };
  
  template<typename t_Object> 
    types::t_unsigned Minimize_Linear<t_Object> :: nb_guess = 0;
  template<typename t_Object> 
    types::t_unsigned Minimize_Linear<t_Object> :: good_guess = 0;
  template<typename t_Object> 
    types::t_unsigned Minimize_Linear<t_Object> :: bad_guess = 0;
  template<typename t_Object> 
    types::t_unsigned Minimize_Linear<t_Object> :: nb_evals = 0;
  template<typename t_Object> 
    types::t_unsigned Minimize_Linear<t_Object> :: nb_grad_evals = 0;

} // namespace opt
#endif // _OPT_MINIMIZE_LINEAR_H_
