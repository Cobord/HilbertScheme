from __future__ import annotations
from typing import Callable, Dict, List, Optional, Set, Tuple
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
        def permute(self,permutation : List[int]):
            if sorted(permutation)==list(range(0,self.my_dimension)):
                self.my_pt = [self.my_pt[idx] for idx in permutation]
            else:
                raise ValueError(f"{permutation} was not a permutation of the {self.my_dimension} coordinate axes")

    AffineDSpacePt.__name__ = 'A^%d (%s)' % (d,field.__name__)
    return AffineDSpacePt

class DPartition(ABC):
    @abstractmethod
    def n_value(self):
        pass
    @abstractmethod
    def d_value(self):
        pass
    @abstractmethod
    def permute(self, permutation : List[int]):
        pass

class TwoDPartition(DPartition):

    @classmethod
    def next_partition(cls,n : int, m : Optional[int] = None) -> TwoDPartition:
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
    def append(self,new_last : int) -> TwoDPartition:
        return TwoDPartition(self.my_partition + [new_last])
    def n_value(self) -> int:
        return self.my_n
    def d_value(self) -> int:
        return 2
    def permute(self,permutation : List[int]):
        if permutation == [0,1]:
            pass
        elif permutation == [1,0]:
            new_partition = [0]*self.my_partition[0]
            for idx in range(0,self.my_partition[0]):
                new_partition[idx] = sum([part > idx for part in self.my_partition])
            self.my_partition = new_partition
        else:
            raise ValueError(f"Only transposing is the proper permutation for a 2D partition")
    def as_list(self) -> List[int]:
        return self.my_partition
    def __str__(self) -> str:
        return str(self.my_partition)

class ThreeDPartition(DPartition):

    @classmethod
    def next_partition(cls,n : int, rst : Optional[Tuple[int,int,int]]) -> ThreeDPartition:
        """Plane Partition of n that fits in box of size (r,s,t)"""
        if n<0:
            raise ValueError(f"n must be non-negative")
        raise NotImplementedError("Generator for all Plane Partitions not implemented")

    def __init__(self,my_partition : Dict[Tuple[int,int],int]):
        all_keys = my_partition.keys()
        row_indices = [tup[0] for tup in all_keys]
        col_indices = [tup[1] for tup in all_keys]
        pi_ij_values = my_partition.values()
        all_natural : Callable[[List[int]],bool] = lambda xs : all([x>=0 for x in xs])
        if not(all_natural(row_indices)) or not(all_natural(col_indices)) or not(all_natural(pi_ij_values)):
            raise ValueError("The i,j for pi_{i,j} must be natural numbers as well as the pi_{i,j} themselves")
        if len(row_indices)==0:
            self.my_partition = [[0]]
            self.my_n = 0
        else:
            max_row = max(row_indices)
            max_col = max(col_indices)
            self.my_partition = [[0]*(max_col+1) for _ in range(max_row+1)]
            for row_idx in range(0,max_row+1):
                for col_idx in range(0,max_col+1):
                    pi_ij = my_partition.get((row_idx,col_idx),0)
                    self.my_partition[row_idx][col_idx] = pi_ij
                    pi_ijp1 = my_partition.get((row_idx,col_idx+1),0)
                    pi_ip1j = my_partition.get((row_idx+1,col_idx),0)
                    if pi_ij<pi_ip1j or pi_ij<pi_ijp1:
                        raise ValueError("pi_{i,j} fails to satisfy the required inequalities that make it interpretable as a plane partition")
            self.my_n = sum(pi_ij_values)
    def n_value(self) -> int:
        return self.my_n
    def d_value(self) -> int:
        return 3
    def __str__(self) -> str:
        return str(self.my_partition)
    def _cyclic_rotate(self) -> None:
        raise NotImplementedError("No cyclic rotation [1,2,0] of plane partition")
    def _xy_swap(self) -> None:
        num_rows = len(self.my_partition)
        num_cols = len(self.my_partition[0])
        new_partition = [[self.my_partition[j][i] for j in range(0,num_rows)] for i in range(0,num_cols)]
        self.my_partition = new_partition
    def permute(self,permutation : List[int]):
        if permutation == [0,1,2]:
            pass
        elif permutation == [1,0,2]:
            self._xy_swap()
        elif permutation == [1,2,0]:
            self._cyclic_rotate()
        elif permutation == [2,0,1]:
            self._cyclic_rotate()
            self._cyclic_rotate()
        elif permutation == [0,2,1]:
            self._cyclic_rotate()
            self._cyclic_rotate()
            self._xy_swap()
        elif permutation == [2,1,0]:
            self._cyclic_rotate()
            self._xy_swap()
        else:
            raise ValueError(f"Only permutations of 0,1,2 are the proper permutation for a 3D partition")

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
        def factory(cls,**kwargs) -> HilbertScheme:
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
        
        def fully_torus_fixed(self) -> bool:
            return all(self.torus_fixed)

        @abstractmethod
        def rotate(self,multiplicative_factors:List[my_field]) -> None:
            pass
        @abstractmethod
        def permute(self,my_permutation:List[int]) -> None:
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
            elif isinstance(partition,ThreeDPartition) and d==3:
                self.partition = partition
            elif partition.d_value() == d:
                raise NotImplementedError(f"d-dimensional partitions for d>3 not supported")
            else:
                raise ValueError(f"Partition's dimension must match with d")

        def __str__(self) -> str:
            return f"The monomial ideal of length {n} for A^{d} ({self.my_field.__name__}) specified by the partition {str(self.partition)}"
        
        def rotate(self,_) -> None:
            return
        def permute(self,permutation : List[int]) -> None:
            if sorted(permutation)==list(range(0,self.my_d)):
                self.partition.permute(permutation)
            else:
                raise ValueError(f"{permutation} was not a permutation of the {self.my_d} coordinate axes")
    
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
                # CAUTION the type of field must have __eq__ such that __eq__(x,0) makes sense
                self.torus_fixed[i] = True if all({pt.my_pt[i]==0 for pt in self.zset}) else False
        
        def __str__(self) -> str:
            return f"The {n} point(s) ({[str(z) for z in self.zset]}) each in {self.point_class.__name__}"
        
        def rotate(self,multiplicative_factors : List[field]) -> None:
            if len(multiplicative_factors) != self.my_d:
                raise ValueError(f"Expected {self.my_d} multiplicative factors")
            # CAUTION the type of field must have __eq__ such that __eq__(x,0) makes sense
            if any([f==0 for f in multiplicative_factors]):
                raise ValueError(f"Multiplicative factors cannot be 0")
            rotation_factor = self.point_class(multiplicative_factors)
            self.zset = {zi*rotation_factor for zi in self.zset}
            return
        def permute(self,permutation : List[int]) -> None:
            if sorted(permutation)==list(range(0,self.my_d)):
                self.zset = {zi.permute(permutation) for zi in self.zset}
            else:
                raise ValueError(f"{permutation} was not a permutation of the {self.my_d} coordinate axes")
            
    return HilbertScheme