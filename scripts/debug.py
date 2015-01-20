from test import setup_script
setup_script(__file__)

"""A place to write rough and ready debug scripts."""

from thompson.automorphism import Automorphism

psi = Automorphism.from_file('psi.aut')
phi = Automorphism.from_file('phi.aut')
print((psi**2).test_conjugate_to(phi) is not None)
psi.test_power_conjugate_to(phi)
