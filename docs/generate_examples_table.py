import os, sys
os.chdir('..')
sys.path.append(os.getcwd())

from thompson import load_example, available_examples, forest, plot

names = sorted(available_examples())
format_string = "{0}\n\tSee the :download:`forest <examples/{0}.pdf>` diagram and function :download:`plot <examples/{0}.svg>`.\n\n\t{1} "

os.chdir('docs')

with open('examples_table.txt', 'wt') as f:
	for name in names:
		aut = load_example(name)
		doc = aut.__doc__.split('\n\n', maxsplit=1)[0].strip()
		base = "examples/" + name
		pdf = base + ".pdf"
		svg = base + ".svg"
		
		if not os.path.exists(pdf):
			forest(aut,
				jobname=base,
				name=r"{{\tiny {}}}".format(name.replace("_", " ")),
				display=False,
				horiz=not name.startswith("periodic_QNB_")
			)
		if not os.path.exists(svg):
			plot(aut, dest=svg,
				diagonal=True,
				display=False
			)
		
		print(format_string.format(name, doc), file=f)