#
#  Version: $Id: _fit.py 1396 2009-11-16 02:59:48Z davezac $
#

""" Implements a secant minimization algorithms.

    Secant bisects an interval until the minimum value of a function is found.
    It works as a somewhat complex iterator. Initialization should happen
    outside the loop:
      1. secant = Secant( interval=(0,10), yvalues=(2, 12), is_integer=True)
      2. for x in secant: 
      3.   ...
      4.   secant.sety(??) # set value at point x.
      5. x, y = secant.result
    Secant will throw if the interval is not monotonic.
    This module has not been debugged.
"""


class Secant:
  def __init__(self, interval, yvalues, is_integer = False, tolerance=1e-8, precision=1e-8):
    """ Initializes a secant module with intervals, functional values for each
        endpoints, and optionally, whether the interval consists of integers only.
    """

    print "Secant module has not been debugged."
    assert interval[0] != interval[1], "Interval has zero length.\n" 

    if self.interval[0] == interval[0]:
      self._first, self._end = (interval[0], yvalue[0]), (interval[1], yvalue[1])
    else:
      self._first, self._end = (interval[1], yvalue[1]), (interval[0], yvalue[0])

    self.result = (None, None)
    self._is_integer = is_integer 
    self._tolerance = tolerance 
    self._precision = precision 
    
    if self.is_integer: # integer means end-points are integers and 1 < tolerance < 2
      self._first[0] = int(self._first[0])
      self._end[0] = int(self._end[0])
      self._precision = 1.5

  def __iter__(self): return self 

  def sety(self, y): self.result[1] = y

  def _check_convergence(self):
    from math import fabs
    if self._end[0] - self._first[0] < self._precision or\
       fabs( self._first[1] - self._end[1] ) < self._tolerance:

      if self._end[0] < self._first[0]: self.result = self._end
      else: self.result = self._first

      raise StopIteration

  def next(self):

    if self.result[0] is None: 
      self._check_convergence() # StopIteration thrown here if at all.
      self.result = ( (self._first[0] + self._end[0])/2, None )
      return self.result[0]

    assert y is not None, "sety has not been called.\n"

    a = self.result[1] < self._end[1] 
    b = self._first[1] < self.result[1] 

    if  a and b:  self._end = self.result[1] # right interval
    elif (not a) and (not b): self.first = self.result[1] # left interval
    else: raise RuntimeError,"Interval is not monotonic.\n" 

    self._check_convergence() # StopIteration thrown here if at all.

    # adjusts  mid term.
    self.mid = ( (self._end[0]-self._first[0])/2, None )

    return self.result[0]


