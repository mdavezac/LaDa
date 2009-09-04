#! /usr/bin/python


def main():

  import sys

  if "-h" in sys.argv or "--help" in sys.argv:
    print ">", sys.argv[0], " filename [default: input]" 
    print "Extracts output from directories listed in filename."
    print "Each line of in filename is a directory containing other directories "\
          "with static vasp calculations."

  filename = "input"
  if len(sys.argv) == 2: filename = str(sys.argv[1])
  elif len(sys.argv) > 2: 
    print "Too many parameters on command-line."
    sys.exit(1)

  file = 0
  try: file = open(filename, "r")
  except IOError: 
    print "Could not open file", filename
    exit(1)

  for line in file: 
    if len(line) == 0: continue
    if re.search( r"^\s*\#", line ): continue;
    for dir in line.split():
      if dir[0] == '#': continue;
      if 
      if not path.exists( 
    extract( line )


if __name__ == "__main__":
  main();
