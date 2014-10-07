from test import setup_script
setup_script(__file__)

from thompson.word import Word

from thompson.examples import example_4_5 as aut



aut.dump_mapping()

basis = aut.domain
print('basis is', basis)

results = {
	Word("x a1 a1 a1", 2, 1) : "left semi-infinite",
	Word("x a1 a2", 2, 1)    : "right semi-infinite",
	Word("x a1 a1 a2", 2, 1) : "complete infinite",
	Word("x a2", 2, 1)       : "incomplete",
}

for word, expected in results.items():
	otype = aut._orbit_type(word, basis)
	print(word, 'should have type', expected)
	print(otype)


# basis = ex.quasinormal_basis()
# print(basis)