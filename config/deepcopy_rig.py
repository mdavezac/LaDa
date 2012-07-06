import quantities
def _rigged_deepcopy(self, memo):
  """ Rig by LaDa so deepcopy of scalars works.
  
      This likely a numpy bug, since deepcopying a scalar array yields a
      builtin type, rather than a scalar array. It is ticket#1176 in numpy bug
      list. 
  """
  from quantities import Quantity
  if len(self.shape) == 0:
    return super(Quantity, self).__deepcopy__(memo) * self.units
  return super(Quantity, self).__deepcopy__(memo)
quantities.Quantity.__deepcopy__ = _rigged_deepcopy

