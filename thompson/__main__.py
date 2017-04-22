#call with python -im thompson
import thompson
from thompson import *
from pprint import pprint

print('The following objects are available from thompson in the current session')
pprint(thompson.__all__)

print("Remember that thompson acts on the right:")
print("    a * b    means     'a, followed by b'")