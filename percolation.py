import numpy as np

def find_threshold_bfs(array):
    first_node = 0
    last_node = len(array) - 1
    probabilities = np.unique(array.ravel())
    low = 0
    high = len(probabilities)

    while high - low > 1:
        i = (high + low) // 2
        prob = probabilities[i]
        copied_array = np.array(array)
        copied_array[copied_array < prob] = 0.0
        if bfs(copied_array, first_node, last_node):
            low = i
        else:
            high = i

    return probabilities[low]


def bfs(graph, source, dest):
    """Perform breadth-first search starting at source. If dest is reached,
    return True, otherwise, return False."""
    # Based on http://www.ics.uci.edu/~eppstein/PADS/BFS.py
    # by D. Eppstein, July 2004.
    visited = set([source])
    nodes = np.arange(0, len(graph))
    stack = [(source, nodes[graph[source] > 0])]
    while stack:
        parent, children = stack[0]
        for child in children:
            if child == dest:
                return True
            if child not in visited:
                visited.add(child)
                stack.append((child, nodes[graph[child] > 0]))
        stack.pop(0)
    return False

