#https://stackoverflow.com/a/13174701/5252017

import sys

def info(type, value, tb):
   if hasattr(sys, 'ps1') or not sys.stderr.isatty():
      # we are in interactive mode or we don't have a tty-like
      # device, so we call the default hook
      sys.__excepthook__(type, value, tb)
   else:
      import traceback, pdb
      # we are NOT in interactive mode, print the exception...
      traceback.print_exception(type, value, tb)
      print
      # ...then start the debugger in post-mortem mode.
      pdb.pm()

sys.excepthook = info

from plmaps   import *
from thompson import *

#1. Load the definition of alpha^q|_D.
#2. Generate lots of elements of the centraliser
#3. Test to see which commute with alpha itself.