from scripts import setup_script, find_examples_passing
setup_script(__file__)

from functools         import lru_cache
from subprocess        import check_output
from random            import randint
from thompson          import *
from thompson.examples import random_automorphism

"""Checks to see if our implementation of Higman's conjugacy test agrees with Hossain's implementation of Belk/Matucci's strand diagrams."""

signature = (2, 1)
def generate_auts():
    """
        >>> phi, psi = generate_auts()
        >>> phi.preserves_order()
        True
        >>> psi.preserves_order()
        True
    """
    phi = random_automorphism(signature, group='F')
    psi = random_automorphism(signature, group='F')
    return phi, psi

def implies(A, B):
    """
        >>> implies(True, True)
        True
        >>> implies(True, False)
        False
        >>> implies(False, True)
        True
        >>> implies(False, False)
        True
    """
    return B or not A

x0 = Automorphism.from_string("""3
(2,1)->(2,1)
x1 a1    -> x1 a1 a1
x1 a2 a1 -> x1 a1 a2
x1 a2 a2 -> x1 a2""")

x1 = Automorphism.from_string("""4
(2,1)->(2,1)
x1 a1       -> x1 a1
x1 a2 a1    -> x1 a2 a1 a1
x1 a2 a2 a1 -> x1 a2 a1 a2
x1 a2 a2 a2 -> x1 a2 a2""")

genrs = [x0, x1]
y0 = ~x0
y1 = ~x1

def test_via_hossain(phi, psi):
    """
        >>> def test(s1, s2): return test_via_hossain(eval_word(s1), eval_word(s2))
        >>> test('x0x0y1x0y1y1y0', 'x0y1x0y1y1')
        True
    """
    phi = as_word(phi) 
    psi = as_word(psi)
    #I think I need to reverse the words for Hossain's tool here
    output = check_output(["java",
      "-jar", r"H:\ConjugacyF\bin\checkConjugacyQuiet.jar",
      phi, psi])
    return output.strip() == b'true'

def as_word(aut):
    """See Cannon, Floyd, Parry 1996: Ex 2.4 and Thm 2.5
    
    >>> f = random_automorphism(signature=(2,1), group='F')
    >>> eval_word(as_word(f)) == f
    True
    >>> for name, aut in {'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1}.items():
    ...     eval_word(name) == aut
    True
    True
    True
    True
    """
    domain_exponents = get_exponents(aut.domain)
    range_exponents  = get_exponents(aut.range)
    # print(domain_exponents)
    # print(range_exponents)
    word = []
    n = len(range_exponents)
    for i in range(n):
        word.append((i,     range_exponents[i]))
    for i in range(n):
        word.append((n-i-1, -domain_exponents[n-i-1]))
    # print(word)
    out = ''.join(expand_letter(*letter) for letter in word)
    return reverse_pair(out)

def reverse_pair(a):
    #http://stackoverflow.com/a/5864313
    return "".join(reversed([a[i:i+2] for i in range(0, len(a), 2)]))

@lru_cache(maxsize=None)
def expand_letter(sub, exp, reversed=False):
    if exp == 0:
        return ""
    letter = 'x' if exp > 0 else 'y'
    if sub in {0, 1}:
        return (letter + str(sub))*abs(exp)
    else:
        #There are two versions of Hossain's applet. The order of composition is switched in them.
        #I don't know what left and right mean in these groups any more.
        # if reversed:
            # return expand_letter(0, sub-1) + expand_letter(1, exp) + expand_letter(0, 1-sub)
        # else:
        return expand_letter(0, 1-sub) + expand_letter(1, exp) + expand_letter(0, sub-1)

def get_exponents(basis):
    return [exponent(leaf, basis) for leaf in basis]

def exponent(leaf, basis):
    exponent = 0
    alphas = [ -x for x in leaf[1:] ]
    while len(alphas) > 1 and alphas.pop() == 1:
        if all(alpha == 2 for alpha in alphas):
            break
        exponent += 1
    return exponent

def eval_word(word):
    """
        >>> eval_word('x0x1x0x1y0x1y1x0') == x0*x1*x0*x1*y0*x1*y1*x0
        True
    """
    basis = Generators.standard_basis((2,1))
    product = Automorphism(basis, basis)
    assert product.is_identity()
    while len(word) >= 2:
        letter, sub, word = word[0], int(word[1]), word[2:]
        aut = genrs[sub]
        if letter == 'y':
            aut = ~aut
        product *= aut
    return product

def search_for_problems():
    num_examples = 0
    num_conjugate_pairs = 0
    while True:
        phi, psi = generate_auts()
        if phi.is_identity() or psi.is_identity():
            continue
        num_examples += 1
        if num_examples % 20 == 0:
            print(num_examples)
        hossain = test_via_hossain(phi, psi)
        me = phi.is_conjugate_to(psi)
        num_conjugate_pairs += me
        if not implies(hossain, me):
            print(phi)
            print(psi)
            phi.save_to_file('phi.aut')
            psi.save_to_file('psi.aut')
            break
    print("tried {} examples. Thompson thought that {} were conjugate".format(num_examples, num_conjugate_pairs))

def debug():
    phi = Automorphism.from_file('phi.aut')
    psi = Automorphism.from_file('psi.aut')
    print(test_via_hossain(phi, psi))
    print(psi.is_conjugate_to(phi))
    globals()['phi'] = phi
    globals()['psi'] = psi

if __name__ == "__main__":
    debug()
    # search_for_problems()
    # import doctest
    # doctest.testmod()
