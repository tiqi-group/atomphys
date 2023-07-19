import numpy as np
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.quantum.cg import CG



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



def is_cyclic_util(graph, v, visited, parent):
    visited[v] = True

    for i in range(len(graph[v])):
        if graph[v, i] != 0:  # if there is a connection
            if visited[i] == False:  # if it has not been visited yet
                if is_cyclic_util(graph, i, visited, v):
                    return True
            elif parent != i:  # if an adjacent vertex is visited and not parent of current vertex
                return True

    return False


def is_cyclic(graph):
    visited = [False] * (len(graph))

    for i in range(len(graph)):
        if visited[i] == False:  # don't recur for v if it is already visited
            if is_cyclic_util(graph, i, visited, -1) == True:
                return True

    return False

import numpy as np
import heapq

def dijkstra(weighted_adj_matrix):
    # Number of nodes
    N = weighted_adj_matrix.shape[0]
    
    # Initialize a list to store the shortest distances
    # Start with all distances as infinity except for the start node
    shortest_distances = [np.inf] * N
    shortest_distances[0] = 0

    # Initialize a priority queue and add the start node
    # Each item in the queue is a tuple of (distance, node)
    # We start with the distance to the start node being 0
    queue = [(0, 0)]

    while queue:
        # Get the node with the smallest distance
        current_distance, current_node = heapq.heappop(queue)

        # If the current distance is larger than the stored distance, this path is not optimal and we can skip it
        if current_distance > shortest_distances[current_node]:
            continue

        # For each neighbor of the current node
        for neighbor in range(N):
            # Skip if there is no edge between current node and neighbor
            if weighted_adj_matrix[current_node, neighbor] == 0:
                continue
            
            # Calculate the distance to this neighbor through the current node
            distance = current_distance + weighted_adj_matrix[current_node, neighbor]

            # If this distance is shorter than the previously known shortest distance to this neighbor
            if distance < shortest_distances[neighbor]:
                # Update the shortest distance to this neighbor
                shortest_distances[neighbor] = distance
                # Add the neighbor to the queue
                heapq.heappush(queue, (distance, neighbor))

    # Replace infinities with -1 to indicate that a node is not reachable
    shortest_distances = [-1 if d == np.inf else d for d in shortest_distances]

    return shortest_distances