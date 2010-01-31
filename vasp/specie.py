""" Defines specie methods. """
class Specie:
  """ Holds atomic specie information:  """
  def __init__(self, symbol, path ):
    """ Initializes a specie.
        @param symbol: is the atomic symbol
        @type symbol: str
        @param path: to the directory with the potcar for this particular atomic types.
          This directory should contain a POTCAR or POTCAR.Z file.
        @type path: str
    """
    import os.path

    self.symbol = symbol
    self.path = os.path.expanduser( path )

