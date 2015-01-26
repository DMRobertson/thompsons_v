from test import setup_script
setup_script(__file__)

"""A place to play around with the package."""

from thompson import *

psi = Automorphism.from_file('psi.aut')
phi = Automorphism.from_file('phi.aut')

print('PSI:', psi)
print('PSI to the -2:', psi ** -2)
print('inv PSI squared:', (~psi) ** 2)

print('Correct answer:', ~psi * ~psi)

print(psi ** -2 == (~psi) ** 2)
print(~psi * ~psi == (~psi) ** 2)
print(~psi == psi ** -1)