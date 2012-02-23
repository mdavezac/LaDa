""" Classes to check pickling in Restart """
from collections import namedtuple
Extract = namedtuple("Extract", ['directory', 'success'])
Vasp = namedtuple("Vasp", ['nonscf'])


