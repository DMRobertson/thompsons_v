"""This was a very hacky extension (read: I didn't know what I was doing) which allows you to refer to \labels declared in a tex source, using whatever the last avialable compilation used as the cross reference name.
"""

import re

from docutils.nodes import Text

refcache = None
reference_finder = re.compile(r"""
	^\\newlabel   #Starts with \newlabel
	\{            #Opening brace
	(             #Start a group
		[^}]+         #As least one non-closing-brace
	)             #End the group
	\}\{\{        #Closing brace, two opening braces
	(             #Start a group
		[^}]+         #As least one non-closing-brace
	)             #End the group
	
""", re.VERBOSE)

def collect_references(aux_file_path):
	global refcache
	refcache = {}
	with open(aux_file_path, 'rt') as f:
		for line in f:
			match = reference_finder.match(line)
			if match is None:
				continue
			label, ref = match.groups()
			refcache[label] = ref

#todo:
	#save the refcache to a pickle object, and load from it if the .aux is not modified
	#expose a config variable for the aux file location

def paperref_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
	try:
		ref = refcache[text]
	except KeyError:
		msg = inliner.reporter.error('Paperref: Could not find a reference to label "%s".' % text, line=lineno)
		prb = inliner.problematic(rawtext, rawtext, msg)
		return [Text(text)], [msg]
	else:
		return [Text(ref)], []

def setup(app):
	try:
		collect_references('paper_references.aux')
	except FileNotFoundError as e:
		app.warn(e)
	else:
		app.add_role('paperref', paperref_role)