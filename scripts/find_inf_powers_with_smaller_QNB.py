from scripts import setup_script
setup_script(__file__)

"""A quick hack searching for counterexamples to AJD's conjecture at the start of sec 6.2."""

from thompson.examples import * 
aut = load_example('power_smaller_QNB')
print('Here is an automorphism and its QNB:')
print(aut)
aut.dump_QNB()
print()

print('Here is its square and the square\'s QNB:')
square = aut * aut
print(square)
square.dump_QNB()

print('Note that the first QNB contains xa211 and xa212, whereas the second contains xa21.')

def search():
    num_examples = 0
    while True:
        p, i = random_automorphism().free_factors()
        if i is None:
            continue
        i = Automorphism(i.domain, i.range)
        
        num_examples += 1
        if True or num_examples % 50 == 0:
            print(num_examples)
        X = i.quasinormal_basis
        
        power = i
        for a in range(10):
            power *= i
            Y = power.quasinormal_basis
            if not all(X.is_above(y) for y in Y):
                print(i)
                print(X)
                print(a+2)
                print(power)
                print(Y)
                raise ValueError()