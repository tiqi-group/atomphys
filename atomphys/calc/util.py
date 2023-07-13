import numpy as np
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.quantum.cg import CG



def spherical_basis_vector(q):
    """
    Function to compute the normalized spherical basis vector c(q).
    
    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding basis vector.
    """
    if q == 1:
        return (-1)*np.array([1, -1j, 0]) / np.sqrt(2)
    elif q == 0:
        return np.array([0, 0, 1])
    elif q == -1:
        return np.array([1, 1j, 0]) / np.sqrt(2)
    else:
        raise ValueError("q should be -1, 0, or 1.")
    
def spherical_basis_second_rank_tensor(q):
    """
    Function to compute the second rank tensor c(q)_ij.

    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding second rank tensor.
    """
    tensor = np.zeros((3, 3), dtype=complex)
    for m1 in range(-1, 2):
        for m2 in range(-1, 2):
            wigner_3j_symbol = w3j(1, 1, 1, m1, m2, -q)
            basis_vector_m1 = spherical_basis_vector(m1)
            basis_vector_m2 = spherical_basis_vector(m2)
            tensor += np.sqrt(10 / 3) * ((-1) ** q) * wigner_3j_symbol * np.outer(basis_vector_m1, basis_vector_m2)
    return tensor
