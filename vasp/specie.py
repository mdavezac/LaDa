class Specie:
  """ Holds atomic specie information: 
       - symbol is the atomic symbol
       - file is the path to the potcar directory.
         it will be assumed to be zipped.
  """

  def __init__(self, symbol, path ):
    """ Initializes a specie.
        _ symbol is the atomic sumbol
        _ path is the path to the potcar directory.
    """
    import os.path

    self.symbol = symbol
    self.path = os.path.expanduser( path )

