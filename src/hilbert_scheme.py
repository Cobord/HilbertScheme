"""
hilbert scheme
cases of
- n disjoint points in A^d
- the torus fixed ideal represented by a d dimensional partition
"""

# type :ignore

from __future__ import annotations
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple
from abc import ABC, abstractmethod

# pylint:disable = invalid-name, superfluous-parens, too-many-statements, undefined-variable

def make_AffineDSpacePt(d: int, field: type = complex):
    """
    return a class whose instances represent
    a point in field^d
    by default it is C^d
    """
    class AffineDSpacePt:  # pylint: disable=redefined-outer-name
        """
        the aforementioned class
        """

        def __init__(self, pt: List[field]) -> None: #type:ignore
            """
            provide the coordinates of the point
            """
            if len(pt) == d and all((isinstance(zi, field) for zi in pt)):
                self.my_pt = pt
                self.my_field = field
                self.my_dimension = d
            else:
                raise ValueError(f"pt must be a list of {d} {field} numbers")

        def __str__(self):
            """
            display the values of all the coordinates
            """
            return str(self.my_pt)

        def __eq__(self, other):
            """
            all the coordinates have to match up and the underlying field and dimension
            """
            return self.my_pt == other.my_pt and\
                self.my_field == other.my_field and\
                    self.my_dimension == other.my_dimension

        def __hash__(self) -> int:
            """
            use the underlying point
            """
            return sum(hash(z) for z in self.my_pt)

        def __add__(self, other):
            """
            coordinatewise addition
            """
            if self.my_field != other.my_field or self.my_dimension != other.my_dimension:
                raise ValueError(
                    "Cannot add points in different fields or dimensions")
            new_point = (self.my_field(zi+other.my_pt[i]) for i, zi in enumerate(self.my_pt))
            return AffineDSpacePt(list(new_point))

        def __mul__(self, other):
            """
            coordinatewise multiplication
            """
            if self.my_field != other.my_field or self.my_dimension != other.my_dimension:
                raise ValueError(
                    "Cannot coordinatewise multiply points in different fields or dimensions")
            new_point = (self.my_field(zi*other.my_pt[i]) for i, zi in enumerate(self.my_pt))
            return AffineDSpacePt(list(new_point))

        def permute(self, permutation: List[int]):
            """
            permute the coordinates
            """

            if sorted(permutation) != list(range(0, self.my_dimension)):
                msg = " ".join([f"{permutation} was not a permutation",
                               "of the {self.my_dimension} coordinate axes"])
                raise ValueError(msg)
            self.my_pt = [self.my_pt[idx] for idx in permutation]

    AffineDSpacePt.__name__ = f"A^{d} ({field.__name__})"
    return AffineDSpacePt


class DPartition(ABC):
    """
    d dimensional partition
    with d=2 being the usual partition,
    d=3 being plane partitions
    """
    @abstractmethod
    def n_value(self):
        """
        the total sum
        """
    @abstractmethod
    def d_value(self):
        """
        how many dimensions
        """
    @abstractmethod
    def permute(self, permutation: List[int]):
        """
        permute the d dimensions,
        if imagine as boxes in an orthant of R^d
        do the corresponding S_d subset O(d)
        and move the boxes accordingly
        """


class TwoDPartition(DPartition):
    """
    Young Tableaux
    """

    @classmethod
    def next_partition(cls, n: int, m: Optional[int] = None):
        """Partition n with a maximum part size of m."""
        if n < 0:
            raise ValueError(f"n must be non-negative not {n}")
        if m is None or m >= n:
            yield TwoDPartition([n])
        top_of_range = n-1 if (m is None or m >= n) else m
        for new_m in range(top_of_range, 0, -1):
            for p in TwoDPartition.next_partition(n-new_m, new_m):
                yield p.append(new_m)

    def __init__(self, my_partition: List[int]):
        """
        provide the list of parts
        """
        self.my_partition = my_partition
        self.my_partition.sort(reverse=True)
        self.my_n = sum(self.my_partition)

    def append(self, new_last: int) -> TwoDPartition:
        """
        append another part
        """
        return TwoDPartition(self.my_partition + [new_last])

    def n_value(self) -> int:
        """
        is a partition of this number
        """
        return self.my_n

    def d_value(self) -> int:
        """
        this is a 2d partition
        """
        return 2

    def permute(self, permutation: List[int]):
        """
        either leave the partition alone for the identity permutation
        or transpose the tableaux for the [1,0] permutation
        """
        if permutation == [0, 1]:
            pass
        elif permutation == [1, 0]:
            new_partition = [0]*self.my_partition[0]
            for idx in range(0, self.my_partition[0]):
                new_partition[idx] = sum(
                    (part > idx for part in self.my_partition))
            self.my_partition = new_partition
        else:
            raise ValueError(
                "Only transposing is the proper permutation for a 2D partition")

    def as_list(self) -> List[int]:
        """
        the list of parts
        """
        return self.my_partition

    def __str__(self) -> str:
        """
        the string representation is the list of parts
        """
        return str(self.my_partition)


class ThreeDPartition(DPartition):
    """
    a plane partition
    """

    @classmethod
    def next_partition(cls, n: int, _rst: Optional[Tuple[int, int, int]]) -> ThreeDPartition:
        """Plane Partition of n that fits in box of size (r,s,t)"""
        if n < 0:
            raise ValueError(f"n must be non-negative not {n}")
        raise NotImplementedError(
            "Generator for all Plane Partitions not implemented")

    def __init__(self, my_partition: Dict[Tuple[int, int], int]):
        """
        initialize with a dictionary whose value for the key ij represents
        the number of boxes above that xy location
        """
        all_keys = my_partition.keys()
        row_indices = [tup[0] for tup in all_keys]
        col_indices = [tup[1] for tup in all_keys]
        pi_ij_values = my_partition.values()
        all_natural: Callable[[Iterable[int]],
                              bool] = lambda xs: all((x >= 0 for x in xs))
        if not (all_natural(row_indices)) or\
                not (all_natural(col_indices)) or\
                not all_natural(pi_ij_values):
            raise ValueError("".join(["The i,j for pi_{i,j} must be natural numbers",
                             "as well as the pi_{i,j} themselves"]))
        if len(row_indices) == 0:
            self.my_partition = [[0]]
            self.my_n = 0
        else:
            max_row = max(row_indices)
            max_col = max(col_indices)
            self.my_partition = [[0]*(max_col+1) for _ in range(max_row+1)]
            for row_idx in range(0, max_row+1):
                for col_idx in range(0, max_col+1):
                    pi_ij = my_partition.get((row_idx, col_idx), 0)
                    self.my_partition[row_idx][col_idx] = pi_ij
                    pi_ijp1 = my_partition.get((row_idx, col_idx+1), 0)
                    pi_ip1j = my_partition.get((row_idx+1, col_idx), 0)
                    if pi_ij < pi_ip1j or pi_ij < pi_ijp1:
                        msg = " ".join(["pi_{i,j} fails to satisfy",
                                       "the required inequalities",
                                         "that make it interpretable as a plane partition"])
                        raise ValueError(msg)
            self.my_n = sum(pi_ij_values)

    def n_value(self) -> int:
        """
        is a plane partition of this number
        """
        return self.my_n

    def d_value(self) -> int:
        """
        is a plane partition
        """
        return 3

    def __str__(self) -> str:
        """
        the string representation is the array of heights
        """
        return str(self.my_partition)

    def _cyclic_rotate(self) -> None:
        """
        cyclically rotate the axes
        should only be accessed via permute
        hence the _
        """
        raise NotImplementedError(
            "No cyclic rotation [1,2,0] of plane partition")

    def _xy_swap(self) -> None:
        """
        swap the x and y axes
        should only be accessed via permute
        hence the _
        """
        num_rows = len(self.my_partition)
        num_cols = len(self.my_partition[0])
        new_partition = [[self.my_partition[j][i]
                          for j in range(0, num_rows)] for i in range(0, num_cols)]
        self.my_partition = new_partition

    def permute(self, permutation: List[int]):
        """
        permute the 3 dimensions,
        """
        if permutation == [0, 1, 2]:
            pass
        elif permutation == [1, 0, 2]:
            self._xy_swap()
        elif permutation == [1, 2, 0]:
            self._cyclic_rotate()
        elif permutation == [2, 0, 1]:
            self._cyclic_rotate()
            self._cyclic_rotate()
        elif permutation == [0, 2, 1]:
            self._cyclic_rotate()
            self._cyclic_rotate()
            self._xy_swap()
        elif permutation == [2, 1, 0]:
            self._cyclic_rotate()
            self._xy_swap()
        else:
            raise ValueError(
                "Only permutations of 0,1,2 are the proper permutation for a 3D partition")


def make_HilbertScheme(n: int, d: int, field: type = complex):
    """
    return a class whose instances specify points in Hilb^n (A^d(field))
    """

    if not (n > 0 and d > 0):
        raise ValueError(f"n (provided : {n}) and d (provided : {d}) must be positive integers")

    class HilbertScheme(ABC):
        """
        Specify a point in Hilb^n (A^d(field)) by either
        a d dimensional partition of n
        or n distinct points each in A^d(field)
        """
        torus_fixed: List[bool] = [False]*d
        my_n: int = n
        my_d: int = d
        my_field: type = field

        @classmethod
        def factory(cls, **kwargs) -> HilbertScheme:
            """
            construct the instance of the appropriate subclass
            depending on whether zset or dpartition is provided
            in the keyword arguments
            """
            if not ((zset := kwargs.get('zset', None)) is None):
                from_set = HilbertSchemeSet(zset)
                if from_set.fully_torus_fixed() and d == 2:
                    return HilbertSchemeTorusFixed(partition=TwoDPartition([1]))
                return from_set
            if not ((dpartition := kwargs.get('dpartition', None)) is None):
                return HilbertSchemeTorusFixed(dpartition)
            msg = " ".join([f"Expected either set of {n} points in F^{d}",
                            "or {d} dimensional partition of {n}"])
            raise ValueError(msg)

        def fully_torus_fixed(self) -> bool:
            """
            is this ideal fixed by the (k^*)^d action
            each of the d entries of self.torus_fixed is for corresponding k^* subset (k^*)^d
            """
            return all(self.torus_fixed)

        @abstractmethod
        def rotate(self, multiplicative_factors: List[my_field]) -> None: #type:ignore
            """
            do the torus action
            """

        @abstractmethod
        def permute(self, my_permutation: List[int]) -> None:
            """
            the induced action from a permutation of d
            """

    class HilbertSchemeTorusFixed(HilbertScheme):
        """Indicate a particular monomial ideal with a d dimensional partition of n"""
        partition : Optional[DPartition] = None

        def __init__(self, partition: DPartition):
            """
            a torus fixed ideal specified by a d dimensional partition
            """
            self.torus_fixed = [True]*d
            if not (partition.n_value() == n):
                raise ValueError(f"Must be a partition of {n}")
            if isinstance(partition, TwoDPartition) and d == 2:
                self.partition = partition
            elif isinstance(partition, ThreeDPartition) and d == 3:
                self.partition = partition
            elif partition.d_value() == d:
                raise NotImplementedError(
                    f"d-dimensional partitions for d>3 not supported, wanted {d}")
            else:
                msg = " ".join([f"Partition's dimension ({partition.d_value()})",
                                "must match with d ({d})"])
                raise ValueError(msg)

        def __str__(self) -> str:
            """
            describe this ideal as a point of the corresponding hilbert scheme
            """
            return " ".join([f"The monomial ideal of length {n} for A^{d}",
                             f"({self.my_field.__name__}) specified by the",
                             f"partition {str(self.partition)}"])

        def rotate(self, _) -> None:
            """
            this is torus fixed so nothing happens
            """
            return

        def permute(self, my_permutation: List[int]) -> None:
            """
            the induced action from a permutation of d
            ends up permuting the DPartition
            """
            if self.partition is None:
                return
            if sorted(my_permutation) == list(range(0, self.my_d)):
                self.partition.permute(my_permutation)
            else:
                raise ValueError(
                    f"{my_permutation} was not a permutation of the {self.my_d} coordinate axes")

    class HilbertSchemeSet(HilbertScheme):
        """
        Indicate a point in the nonsingular locus of the
        Hilbert-Chow morphism by specifying a set of
        n distinct points in A^d(field)
        """

        def __init__(self, my_zset: Set[AffineDSpacePt]): #type:ignore
            """
            given a set of n points in A^d (k)
            get a point in Hilb^n (A^d) (k)
            """
            if len(my_zset) == 0:
                raise ValueError("Expected a list of more than 0")
            expected_field = next(iter(my_zset)).my_field
            if not (self.my_field == expected_field):
                msg = " ".join([f"Wanted the set to be of {self.my_field.__name__}",
                                "numbers rather than {expected_field.__name__}"])
                raise ValueError(msg)
            expected_class = make_AffineDSpacePt(d, field=expected_field)
            if any((not (pt.__class__.__name__ == expected_class.__name__) for pt in my_zset)):
                raise ValueError(
                    f"There was a point not in affine {d} space over {expected_field.__name__}")
            self.point_class = expected_class
            self.zset : Set[expected_class] = set() #type:ignore
            for pt in my_zset:
                self.zset.add(pt) #type:ignore
            if not (len(self.zset) == n): #type:ignore
                raise ValueError(f"Set not of size {n}")
            for i in range(d):
                # CAUTION the type of field must have __eq__ such that __eq__(x,0) makes sense
                self.torus_fixed[i] = all({pt.my_pt[i] == 0 for pt in self.zset}) #type:ignore

        def __str__(self) -> str:
            """
            describe this ideal as a point of the corresponding hilbert scheme
            """
            return " ".join([f"The {n} point(s) ({[str(z) for z in self.zset]})", #type:ignore
                             f"each in {self.point_class.__name__}"])

        def rotate(self, multiplicative_factors: List[field]) -> None: #type:ignore
            """
            the rotation acts on each of the n points
            """
            if len(multiplicative_factors) != self.my_d:
                raise ValueError(
                    f"Expected {self.my_d} multiplicative factors")
            # CAUTION the type of field must have __eq__ such that __eq__(x,0) makes sense
            if any((f == 0 for f in multiplicative_factors)):
                raise ValueError("Multiplicative factors cannot be 0")
            rotation_factor = self.point_class(multiplicative_factors)
            self.zset = {zi*rotation_factor for zi in self.zset} #type:ignore

        def permute(self, my_permutation: List[int]) -> None:
            """
            the induced action from a permutation of d
            ends up permuting each point
            """
            if sorted(my_permutation) == list(range(0, self.my_d)):
                self.zset = {zi.permute(my_permutation) for zi in self.zset}
            else:
                raise ValueError(
                    f"{my_permutation} was not a permutation of the {self.my_d} coordinate axes")

    return HilbertScheme
