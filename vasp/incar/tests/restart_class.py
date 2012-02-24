""" Classes to check pickling in Restart """
from collections import namedtuple
Extract = namedtuple("Extract", ['directory', 'success'])
class Vasp(object):
  def __init__(self, value):
    self.nonscf = value
    self.istart = None
    self.icharg = None


