from scripts import setup_script, choose_from_enum
from enum import Enum
setup_script(__file__)

"""Checks all the examples in thompson. Those which are not revealing wrt the minimal tree pair (or alternatively the tree pair corresp to the quasinormal_basis) are listed and rendered."""

from thompson import *
from thompson.examples import load_all_examples
import tempfile, os, os.path, traceback, shutil, re

class BasisOptions(Enum):
    minimal = 1
    quasinormal = 2

desc = {
    BasisOptions.minimal: 'Use the minimal tree pair describing the automorphism.',
    BasisOptions.quasinormal: 'Use the tree pair corresponding to the minimal expansion of the quasinormal basis.'
}

if __name__ == '__main__':
    print('This script looks through all the examples of thompson for automorphisms with non-revealing tree pairs.')

    basis_option = choose_from_enum(BasisOptions, desc)
    examples = load_all_examples()
    with tempfile.TemporaryDirectory() as dir:
        os.chdir(dir)
        print('\nCreated temporary directory:', dir)
        
        for name in sorted(examples):
            ex = examples[name]
            if basis_option == BasisOptions.minimal:
                domain = None
            elif basis_option == BasisOptions.quasinormal:
                domain = 'wrt QNB'
            result = ex.test_revealing(domain)
            if result is not None:
                print('The automorphism {} of V_{} is not revealing (consider {}).'.format(
                  name, ex.signature, result))
                if input('Render? (y/n) ').strip().lower() == 'y':
                    try:
                        ex.render(name, name=name)
                    except Exception:
                        traceback.print_exc()
                        print('Continuing to next example.')
        input('Press enter to remove temporary files.\nDrawings currently open will not be removed.')
    print('\nTemporary files removed from', dir)
