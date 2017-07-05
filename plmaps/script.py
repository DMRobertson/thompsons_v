#https://stackoverflow.com/a/13174701/5252017

import sys

def info(type, value, tb):
	if hasattr(sys, 'ps1') or not sys.stderr.isatty():
		# we are in interactive mode or we don't have a tty-like
		# device, so we call the default hook
		sys.__excepthook__(type, value, tb)
	else:
		import traceback, pdb
		# we are NOT in interactive mode, print the exception...
		traceback.print_exception(type, value, tb)
		print
		# ...then start the debugger in post-mortem mode.
		pdb.pm()

sys.excepthook = info

from plmaps   import *
from thompson import *
from fractions import Fraction as F
#1. Load the definition of alpha^q|_D.
beta = CPL2(
	[  0   , F(1,4), F(3,8), F(1,2), F(3,4), F( 7,8),   1   ],
	[F(1,2), F(5,8), F(3,4), F(1,1), F(9,8), F(10,8), F(3,2)]
)
print("beta")
print(beta)

beta_squared = beta * beta
print("beta_squared")
print(beta_squared)

bsr = beta_squared.restriction(0, F(1, 2))
print("beta_squared restricted to a fundamental domain for beta (bsr)")
print(bsr)

#2. Generate lots of elements of the centraliser
cent_generator = bsr.one_bump_cent_gen()
print("cent_generator")
print(cent_generator)
assert cent_generator * cent_generator == bsr

#3. Test to see which commute with alpha itself.
#TODO: extend cent_generator to the whole circle
#assert cent_generator.commutes(beta)

assert (~beta * beta).is_identity()

LHS = beta.restriction_of_range(*ends(cent_generator.domain))
prod = LHS * cent_generator * ~beta

prod.commutes(beta_squared.restriction(F(1, 2), 1))

#glue conj and prod together
mapping = { d: r for d, r in cent_generator }
mapping.update( pair for pair in prod )
mapping = sorted(mapping.items())
domain, range = zip(*mapping);
cent_extended = CPL2(domain, range)

print("cent_extended")
print(cent_extended)

assert cent_extended.commutes(beta)
