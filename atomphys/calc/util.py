import numpy as np
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
        return (-1) * np.array([1, -1j, 0]) / np.sqrt(2)
    elif q == 0:
        return np.array([0, 0, 1])
    elif q == -1:
        return np.array([1, 1j, 0]) / np.sqrt(2)
    else:
        raise ValueError(f"Invalid value {q} for q. It must be one of (-1, 0, 1).")


def spherical_basis_second_rank_tensor(q):
    """
    Function to compute the second rank tensor c(q)_ij in cartesian coordinates.

    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding second rank tensor.
    """

    if q == 2:
        return (1 / np.sqrt(6)) * np.array([[1, -1j, 0], [-1j, -1, 0], [0, 0, 0]])
    elif q == 1:
        return (1 / np.sqrt(6)) * np.array([[0, 0, -1], [0, 0, 1j], [-1, 1j, 0]])
    elif q == 0:
        return (1 / 3) * np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 2]])
    elif q == -1:
        return (1 / np.sqrt(6)) * np.array([[0, 0, 1], [0, 0, 1j], [1, 1j, 0]])
    elif q == -2:
        return (1 / np.sqrt(6)) * np.array([[1, 1j, 0], [1j, -1, 0], [0, 0, 0]])
    else:
        raise ValueError(f"Invalid value {q} for q. It must be one of (-2, -1, 0, 1, 2).")


def find_rotating_frame(fields):

    G = nx.Graph()  # create a new graph

    for electric_field, transitions in fields.items():
        for transition in transitions:
            # create an edge between state_i and state_f
            G.add_edge(transition.state_i, transition.state_f,
                       weight=electric_field.angular_frequency.to('MHz').magnitude)

    adjacency_matrix_sparse = nx.adjacency_matrix(G)
    adjacency_matrix = adjacency_matrix_sparse.toarray()
    nodes = list(G.nodes())
    array = adjacency_matrix

    for i in range(array.shape[0]):
        for j in range(array.shape[0]):
            if j > i:
                if array[i, j] != 0:
                    en_i = nodes[i].energy
                    en_j = nodes[j].energy
                    val = array[i, j] * np.sign(en_i - en_j) + array[i, i]
                    if array[j, j] == 0:
                        array[j, j] = val
                    elif val != array[j, j]:
                        raise ValueError("It is impossible to go into RWA Picture")
                    else:
                        array[j, j] = val
    output = []
    for i in range(array.shape[0]):
        output.append(array[i, i])
    return output, nodes
