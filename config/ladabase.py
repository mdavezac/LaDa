try: import lada
except: return
else:
  if "ladabase" not in lada.__all__: return 

OUTCARS_prefix     = 'OUTCARs'
""" Name of the collection of OUTCAR files. """
vasp_database_name = 'cid'
""" Name of the database. """
