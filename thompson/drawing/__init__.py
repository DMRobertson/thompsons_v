from .plot import plot
#from .tpd  import tpd

def display(filepath):
	"""From http://stackoverflow.com/a/435669. Opens the given file with the OS's default application."""
	if sys.platform.startswith('darwin'):
		call(('open', filepath))
	elif os.name == 'nt':
		os.startfile(filepath)
	elif os.name == 'posix':
		call(('xdg-open', filepath))
	else:
		raise NotImplementedError
