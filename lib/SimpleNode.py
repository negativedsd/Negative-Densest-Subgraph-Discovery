class SimpleNode:
    Cdegree = None
    qdegree = None
    # {'name of neighbor':[C*w+ - w-, w+ - q*w-]}
    neighbor_dict = None
    fib_node = None

    def __init__(self):
        self.Cdegree = 0
        self.qdegree = 0
        self.neighbor_dict = {}

    def increase_neighbor(self, name, Cdegree, qdegree):
        if name not in self.neighbor_dict:
            self.neighbor_dict[name] = [Cdegree, qdegree]
        else:
            self.neighbor_dict[name][0] += Cdegree
            self.neighbor_dict[name][1] += qdegree
        self.qdegree += qdegree
        self.Cdegree += Cdegree
