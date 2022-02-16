# DIRECTORY: ~/kpyreb/eSims/MSims/isolateValue.py
#
# Function to isolate the values from a read in file. Returns the value.
# This is used in readIn.py when parsing the starType information.
#
  # Isolates the values:
  # 1. Remove everything after the semicolon
  # 2. Delete whitespace
  # 3. Convert to a float
#  
# ARGUMENT:
#   value - Parsed string to isolate. Value found in convertUnits.py from starTypes dir
#
# Called in readIn.py.
# Author KCT
# HISTORY:
#   30 Jun 2017     Created

def isolateValue(value):
    sep = ';'
    info = float((value.split(sep, 1)[0]).strip())
    return info
