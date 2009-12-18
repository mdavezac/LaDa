class Specie:
  """ Holds atomic specie information: 
      symbol is the atomic sumbol
      file is the path to the potcar directory.
           it will be assumed to be zipped.
  """

  def __init__(self, _symbol, _path ):
    """ Initializes a specie.
        _symbol is the atomic sumbol
        _path is the path to the potcar directory.
    """
    import os.path

    self.symbol = _symbol
    self.path = os.path.expanduser( _path )

