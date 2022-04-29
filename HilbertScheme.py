from typing import List, Optional, Set
from abc import ABC, abstractmethod

def AffineDSpacePt(d:int,field : type = complex):
    class AffineDSpacePt:
        def __init__(self,pt : List[field]) -> None:
            if all([isinstance(zi,field) for zi in pt]) and len(pt)==d:
                self.my_pt = pt
                self.my_field = field
                self.my_dimension = d
            else:
                raise ValueError(f"pt must be a list of {d} {field} numbers")
        def __str__(self):
            return str(self.my_pt)
        def __eq__(self,other):
            return self.my_pt == other.my_pt and self.my_field == other.my_field and self.my_dimension == other.my_dimension
        def __hash__(self) -> int:
            return sum(hash(z) for z in self.my_pt)
        def __add__(self,other):
            if self.my_field == other.my_field and self.my_dimension == other.my_dimension:
                field = self.my_field
                d = self.my_dimension
                return AffineDSpacePt([self.my_field(zi+other.my_pt[i]) for i,zi in enumerate(self.my_pt)])
            else:
                raise ValueError(f"Cannot add points in different fields or dimensions")
        def __mul__(self,other):
            if self.my_field == other.my_field and self.my_dimension == other.my_dimension:
                field = self.my_field
                d = self.my_dimension
                return AffineDSpacePt([self.my_field(zi*other.my_pt[i]) for i,zi in enumerate(self.my_pt)])
            else:
                raise ValueError(f"Cannot coordinatewise multiply points in different fields or dimensions")
    AffineDSpacePt.__name__ = 'A^%d (%s)' % (d,field.__name__)
    return AffineDSpacePt

class DPartition(ABC):
    @abstractmethod
    def n_value(self):
        pass

class TwoDPartition(DPartition):

    @classmethod
    def next_partition(cls,n : int, m : Optional[int] = None) -> 'TwoDPartition':
        """Partition n with a maximum part size of m."""
        if n<0:
            raise ValueError(f"n must be non-negative")
        if m is None or m >= n:
            yield TwoDPartition([n])
        top_of_range = n-1 if (m is None or m >= n) else m
        for new_m in range(top_of_range, 0, -1):
            for p in TwoDPartition.next_partition(n-new_m, new_m):
                yield p.append(new_m)

    def __init__(self,my_partition : List[int]):
        self.my_partition = my_partition
        self.my_partition.sort(reverse=True)
        self.my_n = sum(self.my_partition)
    def append(self,new_last : int) -> 'TwoDPartition':
        return TwoDPartition(self.my_partition + [new_last])
    def n_value(self) -> int:
        return self.my_n
    def __str__(self) -> str:
        return str(self.my_partition)

def HilbertScheme(n:int,d:int,field:type=complex):
    """
    Specify a point in Hilb^n (A^d(field)) by either a d dimensional partition of n or n distinct points each in A^d(field)
    """
    if not(n>0 and d>0):
        raise ValueError(f"n and d must be positive integers")

    class HilbertScheme(ABC):
        torus_fixed : List[bool] = [False]*d
        my_n : int = n
        my_d : int = d
        my_field : type = field
        
        @classmethod
        def factory(cls,**kwargs) -> 'HilbertScheme':
            if not ( (zset := kwargs.get('zset',None)) is None):
                from_set = HilbertSchemeSet(zset)
                if from_set.fully_torus_fixed() and d==2:
                    return HilbertSchemeTorusFixed(dpartition=TwoDPartition([1]))
                else:
                    return from_set
            elif not ( (dpartition := kwargs.get('dpartition',None)) is None):
                return HilbertSchemeTorusFixed(dpartition)
            else:
                raise ValueError(f"Expected either set of {n} points in F^{d} or {d} dimensional partition of {n}")
        
        def fully_torus_fixed(self):
            return all(self.torus_fixed)

        @abstractmethod
        def rotate(self,multiplicative_factors:List[my_field]) -> None:
            pass

    class HilbertSchemeTorusFixed(HilbertScheme):
        """Indicate a particular monomial ideal with a d dimensional partition of n"""
        partition = None
        def __init__(self,partition : DPartition):
            self.torus_fixed = [True]*d
            if not(partition.n_value() == n):
                raise ValueError(f"Must be a partition of {n}")
            if isinstance(partition,TwoDPartition) and d==2:
                self.partition = partition
            else:
                raise NotImplementedError(f"Partition must be a TwoDPartition")

        def __str__(self) -> str:
            return f"The monomial ideal of length {n} for A^{d} ({self.my_field.__name__}) specified by the partition {str(self.partition)}"
        
        def rotate(self,_) -> None:
            return
    
    class HilbertSchemeSet(HilbertScheme):
        """
        Indicate a point in the nonsingular locus of the Hilbert-Chow morphism by specifying a set of n distinct points in A^d(field)
        """
        def __init__(self,my_zset : Set[AffineDSpacePt]):
            if len(my_zset) == 0:
                raise ValueError(f"Expected a list of more than 0")
            else:
                expected_field = next(iter(my_zset)).my_field                
            if not(self.my_field == expected_field):
                raise ValueError(f"Wanted the set to be of {self.my_field.__name__} numbers rather than {expected_field.__name__}")
            expected_class = AffineDSpacePt(d,field=expected_field)
            if any([not(pt.__class__.__name__==expected_class.__name__) for pt in my_zset]):
                raise ValueError(f"There was a point not in affine {d} space over {expected_field.__name__}")
            self.point_class = expected_class
            self.zset = set()
            for pt in my_zset:
                self.zset.add(pt)
            if not(len(self.zset)==n):
                raise ValueError(f"Set not of size {n}")
            for i in range(d):
                self.torus_fixed[i] = True if all({pt.my_pt[i]==0 for pt in self.zset}) else False
        
        def __str__(self) -> str:
            return f"The {n} point(s) ({[str(z) for z in self.zset]}) each in {self.point_class.__name__}"
        
        def rotate(self,multiplicative_factors : List[field]) -> None:
            if len(multiplicative_factors) != self.my_d:
                raise ValueError(f"Expected {self.my_d} multiplicative factors")
            if any([f==0 for f in multiplicative_factors]):
                raise ValueError(f"Multiplicative factors cannot be 0")
            rotation_factor = self.point_class(multiplicative_factors)
            self.zset = {zi*rotation_factor for zi in self.zset}
            return
        
    return HilbertScheme

if __name__ == "__main__":
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
    try:
        _ = z1+w2
    except ValueError as e:
        assert str(e) == "Cannot add points in different fields or dimensions"

    Hilb2C2 = HilbertScheme(2,2)
    p1 = Hilb2C2.factory(zset = {z1,z2})
    print(p1)
    p1.rotate([complex(1,0),complex(0,1)])
    print(f"After rotating by (1,1j) this becomes {p1}")
    print()

    Hilb2C4 = HilbertScheme(2,4)
    p2 = Hilb2C4.factory(zset = {w1,w2,w2copy})
    print(p2)
    assert p2.torus_fixed == [True,False,True,True]
    print()

    Hilb1R2 = HilbertScheme(1,2,float)
    p3 = Hilb1R2.factory(zset={x1,x1,x1})
    print(p3)
    print()

    for cur_partition in TwoDPartition.next_partition(2):
        p4 = Hilb2C2.factory(dpartition=cur_partition)
        print(p4)
    print()

    Hilb5C2 = HilbertScheme(5,2)
    partitions_of_5 = TwoDPartition.next_partition(5)
    for cur_partition in partitions_of_5:
        p5 = Hilb5C2.factory(dpartition=cur_partition)
        print(p5)
    print()

    partitions_of_10 = TwoDPartition.next_partition(10)
    try:
        cur_partition = next(partitions_of_10)
        _ = Hilb5C2.factory(dpartition=cur_partition)
    except ValueError as e:
        assert str(e) == "Must be a partition of 5"

    try:
        _ = HilbertScheme(0,5)
    except ValueError as e:
        assert str(e) == "n and d must be positive integers"