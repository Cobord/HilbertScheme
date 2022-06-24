from hilbert_scheme import TwoDPartition
import numpy as np
from sympy import Matrix
from functools import reduce

Nat = int

def jordan_block(block_size : Nat) -> np.array:
    return np.array([[1 if i==j+1 else 0 for i in range(0,block_size)] for j in range(0,block_size)])

def partition_to_nilmat(partition : TwoDPartition) -> np.array :
    block_sizes = partition.as_list()
    blocks = [jordan_block(block_size) for block_size in block_sizes]
    return reduce(lambda block1,block2 : ksum(block1,block2),blocks)

def ksum(block1 : np.array, block2 : np.array) -> np.array:
    n1,m1 = block1.shape
    n2,m2 = block2.shape
    top_right = np.zeros((n1,m2))
    bottom_left = np.zeros((n2,m1))
    left_side = np.vstack((block1,bottom_left))
    right_side = np.vstack((top_right,block2))
    return np.hstack((left_side,right_side))

def partition_power(partition : TwoDPartition,power : Nat) -> TwoDPartition:
    # J_\lambda^k = P^{-1} J_\mu P
    # where \lambda is partition, power is k
    # and the returned value is \mu
    n = partition.n_value()
    nil_mat = partition_to_nilmat(partition)
    nil_mat = Matrix(nil_mat)
    nil_mat = nil_mat**power
    _, J = nil_mat.jordan_form()
    super_diagonal = [J[i,i+1] for i in range(0,n-1)]
    block_running_sums = [i+1 for i in range(0,n-1) if super_diagonal[i]==0]
    if len(block_running_sums)==0:
        return TwoDPartition([n])
    if (block_running_sums[-1]<n):
        block_running_sums.append(n)
    block_sizes = [block_running_sums[i]-block_running_sums[i-1] for i in range(1,len(block_running_sums))]
    block_sizes.insert(0,block_running_sums[0])
    return TwoDPartition(block_sizes)