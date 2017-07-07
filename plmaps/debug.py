#https://stackoverflow.com/a/13174701/5252017

def info(type, value, tb):
	#if hasattr(sys, 'ps1') or not sys.stderr.isatty():
	#	sys.__excepthook__(type, value, tb)
	#else:
		import traceback, pdb
		traceback.print_exception(type, value, tb)
		pdb.pm()

def debug():
	import sys
	sys.excepthook = info