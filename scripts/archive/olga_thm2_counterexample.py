from thompson import *
f = load_example('not_conjugate_f')
g = load_example('not_conjugate_g')

print('f conj g?', f.is_conjugate_to(g))
h = Automorphism.from_string("""\
2
(2, 1) -> (2, 1)
x a1   -> x a2 
x a2   -> x a1
""")

print('h.order', h.order)
print('f == ~hgh?', f == ~h * g * h)
print('g == hf~h?', g == h * f * ~h)
