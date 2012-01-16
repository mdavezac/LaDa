""" Database related parameters. """
if "ladabase" in globals()["ladamodules"]:
  OUTCARS_prefix     = 'OUTCARs'
  """ Name of the collection of OUTCAR files. """
  vasp_database_name = 'cid'
  """ Name of the database. """
  username  = "Mayeul d'Avezac"
  """ Name with which to tag file in database. """
  # pymongo_username = 'mdadcast'
  # """ Username in the pymongo database. """
  pymongo_host = 'sccdev'
  """ Host of the database. """
  pymongo_port = 27016
  """ Port to which to connect on host. """
  local_push_dir = "/tmp/database_tmp"
  """ Directory where files are pushed, before being pulled to redrock. """
  ladabase_doconnect = False
  """ Whether to connect to database when starting ipython. """
  add_push_magic_function = False
  """ Whether to the %push magic function to the IPython interface. 
  
      Thanks to the paranoids with insecurity issues, it may not be possible to
      access the database from just anywhere. Hence %push is machine dependent.  
  """
