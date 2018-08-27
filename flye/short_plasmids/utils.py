

def find_connected_components(graph):
    def dfs(vertex, connected_components_counter):
        connected_components[vertex] = connected_components_counter
        used[vertex] = True
        for neighbour in graph[vertex]:
            if not used[neighbour]:
                dfs(neighbour, connected_components_counter)

    n_vertices = len(graph)
    connected_components = [0 for _ in range(n_vertices)]
    connected_components_counter = 0
    used = [0 for _ in range(n_vertices)]

    for i in range(n_vertices):
        if not used[i]:
            dfs(i, connected_components_counter)
            connected_components_counter += 1

    return connected_components, connected_components_counter
