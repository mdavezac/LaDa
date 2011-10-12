#! /usr/bin/python
#
#  Version: $Id: sum_of_separables.cc 1329 2009-09-29 22:53:40Z davezac $
#


def leave_many_out( _sumofseps, _collapse, _structures, ratio=0.33333, \
                    nbsets = 3, tolerance=1e-14, itermax=50, verbose = 0 ):

  from math import fabs as abs
  from random import shuffle
  from allsq import allsq
  import numpy

  # a max function we will use later.
  def get_max( _list ):
    result = None
    for i, u in enumerate(_list):
      if result is None: result = (i, u)
      elif u > result[1]: result = (i, u)
    return result

  def errors( _list, _verbose ):
    result = []
    for i in _list: 
       e = _sumofseps( _structures[i] )
       result.append( abs(e-_structures[i].energy) / float(len(_structures[i].atoms)) ) 
    return result

  zero = 1e-10
  full_list = [i for i in range(0, len(_structures)) if _collapse.get_weight(i) > zero]
  nb_out = int(len(full_list)*ratio)
  assert nb_out < len(full_list), "ratio too large."
  if nb_out == 0: nb_out = 1

  print "Starting leave-many-out for %2i of %4i/%4i structures.\n" \
        % (nbsets, nb_out, len(full_list) )

  plefout, allout, allin, vcgs = 1, 2, 3, 4
  if verbose >= plefout: print "Will output results for each left-out set."
  if verbose >= allout: print "Will output results for each structure in each left-out set."
  if verbose >= allin: print "Will output results for each structure in each left-in set."
  if verbose >= vcgs: print "Will make fitting procedure verbose."


  # creates sets of structures which are left out.
  outs = []
  for i in range(nbsets):
    if len(full_list) < nb_out:
      full_list = [i for i in range(0, len(_structures)) if _collapse.get_weight(i) > zero]
      shuffle( full_list )

    outs.append( full_list[:nb_out] )
    full_list = full_list[nb_out:] 

  full_list = [i for i in range(0, len(_structures)) if _collapse.get_weight(i) > zero]
  weights = [_collapse.get_weight(i) for i in full_list]  



  # loop over sets.
  results = []
  for out in outs:

    # changes weights and saves old values.
    reset = [ (u, _collapse.set_weight(u)) for u in out]

    # perform leave fit.
    allsq( _collapse, verbose = verbose >= vcgs, tolerance = tolerance, itermax=itermax )

    # resets weights.
    for u in reset: _collapse.set_weight(u[0], u[1])

    # reassigns coefficients to _sumofseps.
    _collapse.reassign()

    # index of fitting strutures
    in_ = [ i for i in full_list if not (i in out) ]

    # gets errors using local function define above.
    errors_in  = errors( in_, verbose>=allin)
    errors_out = errors( out, verbose>=allout)

    # computes averages
    average_in  = numpy.average( errors_in, weights=[ weights[i] for i in in_ ] )
    average_out = numpy.average( errors_out, weights=[ weights[i] for i in out ] )

    # gets max element.
    max_in  = get_max( weights[i] * errors_in[j] for j, i in enumerate(in_) )
    max_out = get_max( weights[i] * errors_out[j] for j, i in enumerate(out) )

    # appends result
    results.append( (average_in, max_in[1], average_out, max_out[1]) )

    if verbose < plefout: continue

    print "Set: ", out
    print "Average fitting error: %18.12f  Max fitting error %18.12f\n"\
          "Average prediction error: %18.12f  Max prediction error %18.12f\n"\
          % tuple(results[-1])


  print "\nTotal average fitting error: %18.12f  Total max fitting error %18.12f\n"\
        "Average prediction error: %18.12f  Max prediction error %18.12f\n" \
        % (\
            numpy.average([u[0] for u in results]), 
            max( u[1] for u in results), 
            numpy.average([u[2] for u in results]), 
            max( u[3] for u in results)
          )

    








  
