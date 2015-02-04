from pathlib import Path
import re

def bibliography_entries(tex_source):
	with Path(tex_source).open('rt') as f:
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

translator = {
	'``': '"',
	"''": '"',
	"~": ' ',
	r"\'{e}": "é",
	r"\'{o}": "ó",
}

italics = re.compile(r"""
	\{\\it(.+?)\}
""", re.VERBOSE)

def tidy_up(text):
	text = text.strip()
	for key, value in translator.items():
		text = text.replace(key, value)
	m = italics.match(text)
	print(text, m.groups() if m is not None else None)
	text = italics.sub('*\1*', text)
	return text

if __name__ == '__main__':
	refs = {}
	
	for line in bibliography_entries(r'H:\barker\conj_paper\pconj.tex'):
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
			print('.. [{}] {}'.format(key, text), file=f)