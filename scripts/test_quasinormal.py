from test import setup_script
setup_script(__file__)

print('importing example_4_25')
from thompson.examples import example_4_25 as ex
ex.dump_mapping()

print("\nNow for some new computations:")

for word in ["x1 a1 a2", "x1 a1 x1 a1 a2 L", "x1 a1 a1 a2 a1 a2 x1 a2 a2 a2 L"]:
	print("{!s: <37} -> {!s}".format(word, ex.image(word)))

print("\nWe can see that these calculations have been cached if we dump the mapping again.")
ex.dump_mapping()
print("\nHow about the inverse?")
ex.dump_mapping(inverse=True)