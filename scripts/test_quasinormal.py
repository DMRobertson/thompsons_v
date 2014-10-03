from test import setup_script
setup_script(__file__)

from thompson.word import Word

print('importing example_4_25')
from thompson.examples import example_4_25 as ex
ex.dump_mapping()
ex.dump_mapping(inverse=True)

basis = ex.quasinormal_basis()
print(basis)