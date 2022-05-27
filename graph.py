import numpy as np
from geometry import are_collinear, are_coplanar


class PeriodicGraph:

    def __init__(self,
                 vertex_number,
                 edges=[],
                 translations=[],
                 vertex_attributes={},
                 edges_attributes={}):

        self.v_number = vertex_number
        self.edges = np.array(edges)
        self.translations = np.array(translations, dtype=int)
        self.vertex_attributes = {k: np.array(v) for k, v in vertex_attributes.items()}
        self.edges_attributes = {k: np.array(v) for k, v in edges_attributes.items()}
        self.adjacency_list = None
        self.calc_adjacency()
        self.min_path_weight = None
        self.paths = []

    def calc_adjacency(self):

        self.adjacency_list = [[] for _ in range(self.v_number)]
        for i, (v_1, v_2) in enumerate(self.edges):
            self.adjacency_list[v_1].append((v_2, self.translations[i]))
            self.adjacency_list[v_2].append((v_1, -self.translations[i]))
        return self.adjacency_list

    def get_neighbors(self, vertex, translation):
        return [(n, t + translation) for n, t in self.adjacency_list[vertex]]

    def get_vertices(self, selector=lambda v: True):

        return [i for i in range(self.v_number) if selector(i)]

    def get_edges(self, selector=lambda e: True):

        return [i for i in range(self.e_number) if selector(i)]

    def get_edge_induced_subgraph(self, edges_ind):

        edges_ind = list(set(edges_ind))
        translations = self.translations[edges_ind]
        vertices = np.sort(list((set(v for e in self.edges[edges_ind] for v in e))))
        new_indexes = {j: i for i, j in enumerate(vertices)}
        edges = [(new_indexes[i], new_indexes[j]) for i, j in self.edges[edges_ind]]
        edges_attributes = {k: v[edges_ind] for k, v in self.edges_attributes.items()}
        vertex_attributes = {k: v[vertices] for k, v in self.vertex_attributes.items()}
        return PeriodicGraph(len(vertices), edges, translations, vertex_attributes, edges_attributes)

    def get_vertex_induced_subgraph(self, vertices):

        vertices = set(vertices)
        edges_ind = [i for i, (j, k) in enumerate(self.edges) if j in vertices and k in vertices]
        vertices = np.sort(list(vertices))
        translations = self.translations[edges_ind]
        new_indexes = {j: i for i, j in enumerate(vertices)}
        edges = [(new_indexes[i], new_indexes[j]) for i, j in self.edges[edges_ind]]
        edges_attributes = {k: v[edges_ind] for k, v in self.edges_attributes.items()}
        vertex_attributes = {k: v[vertices] for k, v in self.vertex_attributes.items()}
        return PeriodicGraph(len(vertices), edges, translations, vertex_attributes, edges_attributes)

    def calc_periodicity(self):

        traversed = [None for _ in range(self.v_number)]
        for i in range(self.v_number):
            if traversed[i] is None:
                translations = []
                traversed[i] = np.array([0, 0, 0])
                vertices = [(i, np.array([0, 0, 0]))]
                for j, (v, t_1) in enumerate(vertices):
                    #print("fv", v, t_1, self.vertex_attributes["atoms"][v].fract_coord + t_1)
                    for n, t_2 in self.get_neighbors(v, t_1):
                        #print("nv", n, t_2, v, t_1, traversed[n],
                        #      self.vertex_attributes["atoms"][n].fract_coord + t_2,
                        #      self.vertex_attributes["atoms"][v].fract_coord + t_1)
                        translation_image = traversed[n]
                        if translation_image is None:
                            traversed[n] = t_2
                            vertices.append((n, t_2))
                        elif (translation_image != t_2).any():
                            #print(t_2 - translation_image, t_2, translation_image)
                            translations.append(t_2 - translation_image)
                yield vertices, self.get_periodicity(translations)

    @staticmethod
    def get_periodicity(translations):

        translations
        if len(translations) == 0:
            return 0
        elif are_collinear(translations, tol=0):
            return 1
        elif are_coplanar(translations, tol=0):
            return 2
        return 3

    def find_paths(self, n_1, n_2, n_weights, length_limit=100, tol=1e-5):

        def find(path, traversed, weight):

            if sum([1 for t in traversed if t is not None]) < length_limit:
                neighbor_weight = [(n, n_weights[n[0]]) for n in
                                   self.get_neighbors(*path[-1]) if traversed[n[0]] is None]
                neighbor_weight.sort(key=lambda nw: nw[1])
                for n, w in neighbor_weight:
                    if weight < w:
                        weight = w
                    if weight < self.min_path_weight + tol:
                        new_path = path[:]
                        new_path.append(n)
                        if n[0] != n_2[0]:
                            new_traversed = np.array(traversed)
                            new_traversed[n[0]] = n
                            find(new_path, new_traversed, weight)
                        elif all(n[1] == n_2[1]):
                            self.paths.append((new_path, weight))
                            if weight < self.min_path_weight:
                                self.min_path_weight = weight

        n_weights = np.array(n_weights)
        self.paths, self.min_path_weight = [], float("inf")
        traversed = np.full(self.v_number, None)
        traversed[n_1[0]] = n_1
        find([n_1], traversed, n_weights[n_1[0]])
        paths = [(p, w) for p, w in self.paths if w - self.min_path_weight < tol]
        paths.sort(key=lambda x: (x[1], len(x[0])))
        self.paths = [p for p, w in paths]
        return self.paths

    def find_weightless_periodic_subgraphs(self, weights, tol=1e-6):

        subgraphs = [None, None, None]
        edge_weights = [(i, w) for i, w in enumerate(weights)]
        edge_weights.sort(key=lambda x: x[1])
        edges, weights = [i for i, _ in edge_weights], [w for _, w in edge_weights]
        max_weight = weights[-1]
        while len(edges) != 0:
            for w in weights[::-1]:
                if max_weight - w > tol or len(edges) == 0:
                    max_weight = w
                    break
                weights.pop(-1)
                edges.pop(-1)
            if len(edges) != 0:
                max_periodicity = 0
                periodic_subgraph_nodes = []
                subgraph = self.get_edge_induced_subgraph(edges)

                for nodes, periodicity in subgraph.calc_periodicity():
                    if periodicity != 0:
                        if max_periodicity < periodicity:
                            max_periodicity = periodicity
                        periodic_subgraph_nodes.extend([n for n, t in nodes])

                if len(periodic_subgraph_nodes) != 0:
                    p_subgraph = subgraph.get_vertex_induced_subgraph(np.sort(list(set(periodic_subgraph_nodes))))
                    subgraphs[max_periodicity - 1] = p_subgraph
        return subgraphs









