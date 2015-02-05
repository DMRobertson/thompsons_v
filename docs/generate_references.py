from pathlib import Path
import re

def bibliography_entries(tex_source):
	with Path(tex_source).open('rt', encoding='utf-8') as f:
		finished = False
		while not finished:
			line = next(f)
			finished = line.startswith(r'\begin{thebibliography}')
		
		line = next(f)
		finished = False
		while not finished:
			yield line.strip()
			line = next(f)
			finished = line.startswith(r'\end{thebibliography}')

extractor = re.compile(r"""
	\\bibitem	#Bibitem macro
	\[(\w+)\]		#Displayed reference name
	\{(\w+)\}		#Internal label
	(.+)
""", re.VERBOSE)

translator = [
	('``'    , '"'),
	("''"    , '"'),
	("~"     , ' '),
	(r"\'{e}", "é"),
	(r"\'{o}", "ó"),
	("---"   , "—"), #em dash
	("--"    , "–"), #en dash
	(r"\\"    , "")
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
maths = re.compile(r'\$(.+?)\$')

def tidy_up(text):
	text = text.strip()
	for key, value in translator:
		text = text.replace(key, value)
	
	for italic in italics:
		text = italic.sub(r' *\1* ', text)
	for bold in bolds:
		text = bold.sub(r' **\1** ', text)
	text = maths.sub(r' :math:`\1` ', text)
	return text

if __name__ == '__main__':
	refs = {}
	
	for line in bibliography_entries(r'..\..\barker\conj_paper\pconj.tex'):
		if line == '' or line.startswith('%'):
			continue
		try:
			label, _, text = extractor.match(line).groups()
			refs[label] = tidy_up(text)
		except AttributeError:
			print('problem with:', repr(line))
			raise
	from pprint import pprint
	# pprint(refs)
	
	with open('paper_bibliography.txt', 'wt', encoding='utf-8') as f:
		for key, text in sorted(refs.items()):
			print(r'.. [{}] \{}'.format(key, text), file=f)

