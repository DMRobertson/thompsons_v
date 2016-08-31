from pathlib import Path
from os.path import expanduser
import re


def bibliography_entries(tex_source):
	tex_source = expanduser(tex_source)
	output = ""
	with Path(tex_source).open('rt', encoding='utf-8') as f:
		at_bib_start = False
		while not at_bib_start:
			line = next(f)
			at_bib_start = line.startswith(r'\begin{thebibliography}')
		
		line = next(f)
		at_bib_end = False
		while not at_bib_end:
			output += line
			line = next(f)
			at_bib_end = line.startswith(r'\end{thebibliography}')
	return output

extractor = re.compile(r"""
	\\bibitem             #Bibitem macro
	%                     #Commented out
	\[(\w+)\]             #The old displayed reference name
	\n                    #New line
	\{(\w+)\}             #Internal label
	\s?                   #Optional space
	(                     #The citation text
		(?:               #Consists of the following group repeated at least once
			(?!\\bibitem) #Does not begin with another bibitem macro
			.*\n          #Rest of the current line
		)+
	)               
""", re.VERBOSE)

translator = [
	('``'    , '"'),
	("''"    , '"'),
	("~"     , ' '),
	(r"\'{e}", "é"),
	(r"\'{o}", "ó"),
	(r"\"{o}", "ö"),
	(r"{\v\i}" , "i" ),
	("---"   , "—"), #em dash
	("--"    , "–"), #en dash
	(r"\\"   , "" ),
]

tex_group = r"\{{\\{}\s*(.+?)\s*\}}"
tex_macro = r"\\{}\{{(.+?)\}}"
italics = [
	re.compile(tex_group.format('it')),
	re.compile(tex_group.format('em')),
	re.compile(tex_macro.format('emph'))
]
bolds = [
	re.compile(tex_group.format('bf')),
	re.compile(tex_macro.format('textbf'))
]
code = re.compile(tex_group.format('tt'))
maths = re.compile(r'\$(.+?)\$')
url = re.compile(tex_macro.format('url'))

def tidy_up(text):
	text = text.strip()
	if text.startswith('%'):
		return ''
	for key, value in translator:
		text = text.replace(key, value)
	
	for italic in italics:
		text = italic.sub(r' *\1* ', text)
	for bold in bolds:
		text = bold.sub(r' **\1** ', text)
	text = code.sub(r'``\1``', text)
	text = url.sub(r'\1', text)
	text = maths.sub(r' :math:`\1`', text)
	return text

if __name__ == '__main__':
	refs = {}
	bibliography = bibliography_entries('~/barker/conj_paper/pconj_ijac.tex')
	citations = extractor.findall(bibliography)
	
	for ext_name, int_name, text in citations:
		lines = [tidy_up(line) for line in text.split('\n')]
		refs[ext_name] = " ".join(lines)
	
	with open('paper_bibliography.txt', 'wt', encoding='utf-8') as f:
		for key, text in sorted(refs.items()):
			print(r'.. [{}] \{}'.format(key, text), file=f)

