"""
test word functions dealing with crystal combinatorics
"""

import unittest
import random

# pylint:disable = wrong-import-position,import-error,unused-variable, invalid-name, missing-function-docstring,line-too-long,missing-class-docstring

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")
from word_funcs import is_yamanouchi, crystal_e, crystal_f

class TestWordFunc(unittest.TestCase):
    def test_highest_weight_repeated(self):
        for _ in range(100):
            self.yamanouchi_e_defined()

    def yamanouchi_e_defined(self):
        max_n = 4
        a_word = list(random.choices(range(0, max_n+1), k=random.randint(1,10)))
        if is_yamanouchi(a_word):
            for i in range(max_n):
                ei_a_word = crystal_e(a_word,i)
                self.assertEqual(ei_a_word,None, f"{a_word} is Yamanouchi so it should have been highest weight and E_{i} undefined on it, but got {ei_a_word} instead")
        else:
            found_defined_ei = False
            for i in range(max_n):
                found_defined_ei = crystal_e(a_word,i) is not None
                if found_defined_ei:
                    break
            self.assertTrue(found_defined_ei,f"{a_word} is not Yamanouchi so it should have had at least one E_{i} defined on it")

    def round_trip(self):
        max_n = 4
        a_word = list(random.choices(range(0, max_n+1), k=random.randint(1,10)))
        for i in range(max_n):
            fi_a_word = crystal_f(a_word,i)
            if fi_a_word is not None:
                eifi_a_word = crystal_e(fi_a_word,i)
                self.assertEqual(eifi_a_word,a_word,"ei(fi(y)) = y on all y where fi is defined")
        a_word = list(random.choices(range(0, max_n+1), k=random.randint(1,10)))
        for i in range(max_n):
            ei_a_word = crystal_e(a_word,i)
            if ei_a_word is not None:
                fiei_a_word = crystal_f(ei_a_word,i)
                self.assertEqual(fiei_a_word,a_word,"fi(ei(x)) = x on all x where ei is defined")

    def test_round_trip_repeated(self):
        for _ in range(100):
            self.round_trip()

if __name__ == '__main__':
    unittest.main()
