import unittest

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")
from hilbert_scheme import AffineDSpacePt,TwoDPartition,HilbertScheme

class TestHilbert(unittest.TestCase):

    def test_affine(self):
        Aff2 = AffineDSpacePt(2)
        z1 = Aff2([complex(3,4),complex(5,6)])
        Aff4 = AffineDSpacePt(4)
        w1 = Aff4([complex(0,0),complex(1,0),complex(0,0),complex(0,0)])
        z2 = Aff2([complex(2,-1),complex(-4,1)])
        Aff2R = AffineDSpacePt(2,float)
        x1 = Aff2R([5.0,6.5])
        w2 = Aff4([complex(0,0),complex(0,1),complex(0,0),complex(0,0)])
        w2copy = Aff4([complex(0,0),complex(0,1),complex(0,0),complex(0,0)])

        _ = z1+z2
        _ = w1+w2
        with self.assertRaises(ValueError) as cm:
            _ = z1 + w2
            the_exception = cm.exception
            self.assertEqual(str(the_exception),"Cannot add points in different fields or dimensions", "There should have been an error when trying to add points in different affine spaces")

    def test_H2C2_set(self):
        Aff2 = AffineDSpacePt(2)
        z1 = Aff2([complex(3,4),complex(5,6)])
        z2 = Aff2([complex(2,-1),complex(-4,1)])
        Hilb2C2 = HilbertScheme(2,2)
        p1 = Hilb2C2.factory(zset = {z1,z2})
        self.assertEqual(str(p1),"The 2 point(s) (['[(2-1j), (-4+1j)]', '[(3+4j), (5+6j)]']) each in A^2 (complex)")
        p1.rotate([complex(1,0),complex(0,1)])
        self.assertEqual(str(p1),"The 2 point(s) (['[(3+4j), (-6+5j)]', '[(2-1j), (-1-4j)]']) each in A^2 (complex)")
    
    def test_H2C2_partition(self):
        Hilb2C2 = HilbertScheme(2,2)
        for cur_partition in TwoDPartition.next_partition(2):
            p4 = Hilb2C2.factory(dpartition=cur_partition)

    def test_H2C4_set(self):
        Hilb2C4 = HilbertScheme(2,4)
        Aff4 = AffineDSpacePt(4)
        w1 = Aff4([complex(0,0),complex(1,0),complex(0,0),complex(0,0)])
        w2 = Aff4([complex(0,0),complex(0,1),complex(0,0),complex(0,0)])
        w2copy = Aff4([complex(0,0),complex(0,1),complex(0,0),complex(0,0)])
        p2 = Hilb2C4.factory(zset = {w1,w2,w2copy})
        self.assertEqual(p2.torus_fixed,[True,False,True,True],"Rotations of the 1st,3rd,4th coordinates fix this point")
    
    def test_H1R2_set(self):
        Hilb1R2 = HilbertScheme(1,2,float)
        Aff2R = AffineDSpacePt(2,float)
        x1 = Aff2R([5.0,6.5])
        p3 = Hilb1R2.factory(zset={x1,x1,x1})
    
    def test_H5C2_partition(self):
        Hilb5C2 = HilbertScheme(5,2)
        partitions_of_5 = TwoDPartition.next_partition(5)
        for cur_partition in partitions_of_5:
            p5 = Hilb5C2.factory(dpartition=cur_partition)
    
    def test_wrong_n_partition(self):
        Hilb5C2 = HilbertScheme(5,2)
        partitions_of_10 = TwoDPartition.next_partition(10)
        with self.assertRaises(ValueError) as cm:
            cur_partition = next(partitions_of_10)
            _ = Hilb5C2.factory(dpartition=cur_partition)
            the_exception = cm.exception
            self.assertEqual(str(e),"Must be a partition of 5", "There should have been an error when trying to use a partition of 10 for Hilb^5 (A^2 (C))")
    
    def test_zero_edge_case(self):
        with self.assertRaises(ValueError) as cm:
            _ = HilbertScheme(0,5)
            the_exception = cm.exception
            self.assertEqual(str(e),"n and d must be positive integers")

if __name__ == '__main__':
    unittest.main()