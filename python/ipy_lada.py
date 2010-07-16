""" IPython functions and data. """
from .. import jobs

current = None
""" Current dictionary. """

pickle_filename = None
""" Current pickle filename. """

pickle_directory = None
""" Directory of current pickle. """

def explore(self, arg):
  print arg
