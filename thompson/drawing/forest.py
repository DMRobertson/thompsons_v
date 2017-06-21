from functools import partial

from ..generators import Generators
from ..utilities  import intersection_of_trees, intersection_from_domain, handle_domain

#First time setup
template = None
def setup():
	from jinja2 import Environment, PackageLoader
	global template
	if template is not None:
		return template
	env = Environment(
		loader=PackageLoader('thompson', 'drawing'),
		lstrip_blocks = True,
		trim_blocks   = True,
	)
	template = env.get_template('forest_diagram.tpl')
	return template

def forest_code(aut,
	name='',
	domain='wrt QNB',
	include_styles = True,
	horiz=True,
	LTR=True,
	standalone=True,
	draw_revealing=True,
):
	"""
	Generate TikZ code for representing an automorphism *aut* as a forest pair diagram. The code is written to *dest*, which must be a writeable file-like object in text mode.
	
	:param str name: The label used for the arrow between domain and range forests. This is passed directly to TeX, so you can include mathematics by delimiting it with dollars. Note that backslashes are treated specially in Python unless you use a *raw string*, which is preceeded with an ``r``. For instance, try ``name=r'\gamma_1``
	:param `~thompson.generators.Generators` domain: By default, we use the :meth:`minimal expansion <thompson.generators.Generators.minimal_expansion_for>` of the :meth:`quasi-normal basis <thompson.automorphism.Automorphism.compute_quasinormal_basis>` as the leaves of the domain forest. This can be overridden by providing a *domain* argument.
	:param bool include_styles: Should the styling commands be added to this document?
	:param bool horiz: If True, place the range forest to the right of the domain. Otherwise, place it below.
	:param bool LTR: if True and if horiz is True, draw the domain tree to the left of the range tree. If horiz is True but LTR is false, draw the domain tree to the right of the range tree.
	:param bool standalone: If True, create a standalone LaTeX file. Otherwise just create TikZ code.
	:param bool draw_revealing: Should attractor/repeller paths be highlighted in red? 
	"""
	if name.startswith('$') and not name.endswith('$'):
		raise ValueError("Arrow names must end with a $ if they begin with a $.")
	
	#1. Decide which domain to use for plotting.
	domain = handle_domain(domain, aut)
	range, intersection = intersection_from_domain(domain, aut)
	is_repatt_specialised = partial(is_repatt, intersection=intersection, aut=aut)
	
	# Number the leaves.
	domain = [  ( w, i, draw_revealing and is_repatt_specialised(w, True) )
		for (i, w) in enumerate(domain, start=1)]
	range  = [  ( w, i, draw_revealing and is_repatt_specialised(w, False))
		for (i, w) in enumerate(range, start=1) ]
	#Order the range using Higman's words
	range.sort()
	
	template = setup()
	return template.render(
		#options
		name = name,
		domain = domain,
		range = range,
		horiz = horiz,
		standalone = standalone,
		include_styles = include_styles,
		write_word = partial(write_word, intersection = intersection),
		LTR = LTR
	)

def write_word(word, index, highlight, intersection):
	below_intersection = False
	lines = ["root"]
	output = ""
	if len(word) == 1:
		lines[0] += " [point, label=below:\\strut$1$]"
	for subword in word.subwords(discard_root=True):
		lines.append( write_subword(
			subword, index, highlight, below_intersection, subword == word
		))
		if not below_intersection:
			below_intersection = subword in intersection
	return "\n\t\t\t\t".join(lines)

def write_subword(subword, index, highlight, below_intersection, last):
	options = {}
	output = "-- "
	if below_intersection:
		output += "[component] "
		if highlight:
			#TODO I'm using two different ways to style an edge here---doesn't make sense.
			options["target edge style"] = "spine"
	if last:
		options["label"] = "below:\\strut${}$".format(index)
		if highlight:
			options["label"] = "{[repatt label]" + options["label"] + "}"
			options["repatt"] = None
	output += name(subword, options)
	if last:
		output += ","
	return output

def name(word, options):
	"""Introduce a node with the given key/value options."""
	label = word.address()
	pairs = []
	for key, value in options.items():
		if value is None:
			pairs.append(key)
		else:
			pairs.append('{}={}'.format(key, value))
	if pairs:
		label += " [{}]".format( ', '.join(pairs) )
	return label

def is_repatt(word, for_domain, intersection, aut):
	if word in intersection:
		return False
	ctype = aut.orbit_type(word, intersection)[0]
	if not ctype.is_type_B():
		return False
	if (ctype.characteristic[0] < 0) != for_domain:
		return False
	#TODO check that the root you stop at is above where you started
	return True
