from test import setup_script
setup_script(__file__)

"""A while ago I was considering a question in http://arxiv.org/abs/0911.0979 by Bleak and Salazar-Diaz (Question 1). Would the general Python script you have created to work with V be able to look at the orbits of the two elements \psi_1 and \psi_3 in the note attached (written very badly), and their products?"""

from thompson import *

psi_1 = Automorphism.from_file('nathan/psi_1.aut')
psi_3 = Automorphism.from_file('nathan/psi_3.aut')

print('\nPSI_1')
print(psi_1)
psi_1.dump_QNB()

print('\nPSI_3')
print(psi_3)
psi_3.dump_QNB()

product = psi_1 * psi_3
print('\nproduct = psi_1 * psi_3')
print(product)
product.dump_QNB()

product_2 = psi_3 * psi_1
print('\nproduct_2 = psi_3 * psi_1')
print(product_2)
product.dump_QNB()