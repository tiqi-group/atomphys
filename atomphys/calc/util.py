import numpy as np
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.quantum.cg import CG
import networkx as nx



def spherical_basis_vector(q):
    """
    Function to compute the normalized spherical basis vector c(q) in cartesian corrindates..
    
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
    Function to compute the second rank tensor c(q)_ij in cartesian coordinates.

    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding second rank tensor.
    """
    tensor = np.zeros((3, 3), dtype=complex)
    for m1 in range(-1, 2):
        for m2 in range(-1, 2):
            wigner_3j_symbol = float(w3j(1, 1, 1, m1, m2, -q))
            basis_vector_m1 = spherical_basis_vector(m1)
            basis_vector_m2 = spherical_basis_vector(m2)
            tensor += np.sqrt(10 / 3) * ((-1) ** q) * wigner_3j_symbol * np.outer(basis_vector_m1, basis_vector_m2)
    return tensor




def find_rotating_frame(fields):

    G = nx.Graph() # create a new graph

    for electric_field, transitions in fields.items():
        for transition in transitions:
        # create an edge between state_i and state_f
            G.add_edge(transition.state_i, transition.state_f, weight=electric_field.angular_frequency.to('MHz').magnitude)


    adjacency_matrix_sparse = nx.adjacency_matrix(G)
    adjacency_matrix = adjacency_matrix_sparse.toarray()
    nodes = list(G.nodes())
    print(nodes)
    array = adjacency_matrix
        
    for i in range(array.shape[0]):
        for j in range(array.shape[0]):
            if j>i:
                if array[i,j]!=0:
                    en_i = nodes[i].energy
                    en_j = nodes[j].energy
                    val = array[i,j]*np.sign(en_i-en_j)+array[i,i]
                    if array[j,j]==0:
                        array[j,j]=val
                    elif val!=array[j,j]:
                        raise ValueError("It is impossible to go into RWA Picture")
                    else:
                        array[j,j]=val
    output = []
    for i in range(array.shape[0]):
        output.append(array[i,i])
    return output, nodes