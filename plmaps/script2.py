from thompson import *

def dump(aut, name):
	print(name, aut.tikz_path(), sep="\n")

alpha = Automorphism.from_string("""6
(2, 1) -> (2, 1)
x1 a1 a1    -> x1 a1 a1 a1
x1 a1 a2 a1 -> x1 a1 a1 a2
x1 a1 a2 a2 -> x1 a1 a2   
x1 a2 a1    -> x1 a2 a1 a1
x1 a2 a2 a1 -> x1 a2 a1 a2
x1 a2 a2 a2 -> x1 a2 a2   """)

dump(alpha, "alpha")

rot = Automorphism.rotation(1, 2)
dump(rot, "rot")

print("alpha commutes with rot:", alpha.commutes(rot))

beta = alpha * rot
dump(beta, "beta")

beta_squared = beta * beta
dump(beta_squared, "beta_squared")
