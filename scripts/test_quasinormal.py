from test import setup_script
setup_script(__file__)

from thompson import *
from thompson.examples import *
from thompson.orbits import *

psi = example_5_12_psi
phi = example_5_12_phi

psi_p = example_5_12_psi_p
phi_p = example_5_12_phi_p

print('*PSI DIRECTLY PERIODIC*', psi_p)
print()
print('*PHI DIRECTLY PERIODIC*', phi_p)
print()

print('Is psi_p conjugate to phi_p?\n')
rho_p = psi_p.test_conjugate_to(phi_p)
print('*RHO DIRECTLY PERIODIC*', rho_p)
print()

print('Is psi conjugat to phi?\n')
rho = psi.test_conjugate_to(phi)
print('rho:', rho)
print()