from random import random


class SyntheticNode:
    degree = None
    total_degree = None
    neighbor_dict = None

    def __init__(self, n):
        self.degree = [0] * n
        self.neighbor_dict = {}
        self.total_degree = 0

    #     type is int from 0 to len(degree)-1
    def increase_neighbor(self, name, type, degree):
        if name not in self.neighbor_dict:
            self.neighbor_dict[name] = {type: degree}
        else:
            if type not in self.neighbor_dict[name]:
                self.neighbor_dict[name][type] = degree
                self.degree[type] += degree
            else:
                self.neighbor_dict[name][type] = degree


def erdos_renyi(n, p, ratio):
    node_dict = {}
    for i in range(n):
        node_dict[i] = SyntheticNode(2)
    for i in range(n):
        for j in range(i+1, n):
            if random() < p*ratio:
                node_dict[i].increase_neighbor(j, 0, 1)
                node_dict[j].increase_neighbor(i, 0, 1)

    for i in range(n):
        for j in range(i+1, n):
            if random() < p:
                node_dict[i].increase_neighbor(j, 1, -1)
                node_dict[j].increase_neighbor(i, 1, -1)
            # neg = random()*2*p
            # # neg = 1
            # node_dict[i].increase_neighbor(j, 1, -neg)
            # node_dict[j].increase_neighbor(i, 1, -neg)
    return node_dict


def plant(node_dict, n1, p1, n2, p2, clique = True):
    if n1+n2 > len(node_dict):
        return False

    for i in range(n1):
        for j in range(1,p1+1):
            tempIndex = (i+j)%n1
            node_dict[i].increase_neighbor(tempIndex, 0, 1)
            node_dict[tempIndex].increase_neighbor(i, 0, 1)
        # for j in range(i+1, n1):
        #     if clique:
        #         node_dict[i].increase_neighbor(j, 0, 1)
        #         node_dict[j].increase_neighbor(i, 0, 1)
        #     else:
        #         if random() < p1:
        #             node_dict[i].increase_neighbor(j, 0, 1)
        #             node_dict[j].increase_neighbor(i, 0, 1)
    for i in range(len(node_dict) - n2, len(node_dict)):
        for j in range(1,p2+1):
            tempIndex = (len(node_dict) - n2) + (i+j)%n2
            node_dict[i].increase_neighbor(tempIndex, 0, 1)
            node_dict[tempIndex].increase_neighbor(i, 0, 1)
        # for j in range(i+1, len(node_dict)):
        #     if clique:
        #         node_dict[i].increase_neighbor(j, 0, 1)
        #         node_dict[j].increase_neighbor(i, 0, 1)
        #     else:
        #         if random() < p2:
        #             node_dict[i].increase_neighbor(j, 0, 1)
        #             node_dict[j].increase_neighbor(i, 0, 1)

    return node_dict

# test
