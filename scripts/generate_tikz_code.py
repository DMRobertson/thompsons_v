"""I wrote this as part of other work. I thought it should stay as an example!"""

from thompson          import Automorphism
from thompson.examples import olga_f as aut

PATH = r'H:\robertson\Reading\Salazar-Díaz revealing pairs\olga_f_{}.tikz'
A = aut.domain.copy()
aut.write_tikz_code(PATH.format('AB'), A, name='f')

C = A.copy().expand(-1).expand(2)
aut.write_tikz_code(PATH.format('CD'), C, name='f')

E = A.copy().expand(-4).expand(0)
aut.write_tikz_code(PATH.format('EF'), E, name='f')

G = A.copy().expand(0).expand(1)
aut.write_tikz_code(PATH.format('GH'), G, name='f')