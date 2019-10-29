class SimpleNode:
    C_degree = None
    positive_degree = None
    # {'name of neighbor':[C*w+ - w-, w+]}
    neighbor_dict = None
    fib_node = None

    def __init__(self):
        self.C_degree = 0
        self.positive_degree = 0
        self.neighbor_dict = {}

    def increase_neighbor(self, name, Cdegree, positive_degree):
        if name not in self.neighbor_dict:
            self.neighbor_dict[name] = [Cdegree, positive_degree]
        else:
            self.neighbor_dict[name][0] += Cdegree
            self.neighbor_dict[name][1] += positive_degree
        self.positive_degree += positive_degree
        self.C_degree += Cdegree
