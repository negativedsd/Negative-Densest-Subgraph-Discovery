from lib.fibheap import FibonacciHeap
from lib.SimpleNode import SimpleNode
import json


def peeling(node_dict, total_C_degree, total_positive_degree, fib_heap, q, B, C, lambda1, lambda2):
    n = node_dict.__len__()
    C_average_degree = float(total_C_degree) / n
    positive_avg_degree = total_positive_degree / n
    S_size = n
    subgraph = None
    for i in range(n - 1):

        if i % 50000 == 0:
            print(i, C_average_degree, positive_avg_degree, (C * positive_avg_degree - C_average_degree) / (B * q))

        # find min node from graph (remove from heap)
        node_to_remove = fib_heap.extract_min().value
        for neighbor in node_dict[node_to_remove].neighbor_dict.keys():

            # get dictionary that has all edges between two nodes
            C_degree_loss = node_dict[node_to_remove].neighbor_dict[neighbor][0]
            node_dict[neighbor].C_degree -= C_degree_loss
            pos_degree_loss = node_dict[node_to_remove].neighbor_dict[neighbor][1]
            node_dict[neighbor].positive_degree -= pos_degree_loss

            # here the key can be actually increased
            if neighbor != node_to_remove:
                fib_heap.decrease_key(node_dict[neighbor].fib_node, node_dict[neighbor].C_degree)
                del node_dict[neighbor].neighbor_dict[node_to_remove]
            total_C_degree -= C_degree_loss
            total_positive_degree -= pos_degree_loss

        del node_dict[node_to_remove]
        C_average_degree = float(total_C_degree) / (n - i - 1)
        positive_avg_degree = total_positive_degree / (n - i - 1)
        S_size = n - i - 1
        if C_average_degree - (C - 1) * positive_avg_degree > q * lambda2 - lambda1:
            if len(node_dict) < 30:
                subgraph = list(node_dict)
            return True, C_average_degree, positive_avg_degree, S_size, subgraph

    return False, C_average_degree, positive_avg_degree, S_size, subgraph


class RiskNode:
    degree = None
    total_degree = None
    neighbor_dict = None
    paper_count = None

    def __init__(self, n):
        self.degree = [0] * n
        self.neighbor_dict = {}
        self.total_degree = 0
        self.paper_count = 0

    #     type is int from 0 to len(degree)-1
    def increase_neighbor(self, name, type, degree):
        if name not in self.neighbor_dict:
            self.neighbor_dict[name] = {type: degree}
        else:
            if type not in self.neighbor_dict[name]:
                self.neighbor_dict[name][type] = degree
            else:
                self.neighbor_dict[name][type] += degree
        self.degree[type] += degree

    def set_neighbor_risk(self, name, degree):
        self.neighbor_dict[name][1] = degree
        self.degree[1] += degree


def process_tmdb_file(file_path):
    actor_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:
        # print(len(relation_list))
        weight = relation['popularity'] * relation['possibility']
        risk = relation['popularity'] * relation['possibility'] * (1 - relation['possibility'])
        if risk == 0:
            continue
        if relation['actors'][0] not in actor_dict:
            actor_dict[relation['actors'][0]] = RiskNode(2)
        if relation['actors'][1] not in actor_dict:
            actor_dict[relation['actors'][1]] = RiskNode(2)
        actor_dict[relation['actors'][0]].increase_neighbor(relation['actors'][1], 0, weight)
        actor_dict[relation['actors'][1]].increase_neighbor(relation['actors'][0], 0, weight)
        actor_dict[relation['actors'][0]].set_neighbor_risk(relation['actors'][1], -risk)
        actor_dict[relation['actors'][1]].set_neighbor_risk(relation['actors'][0], -risk)
    return actor_dict


def process_dblp_file(file_path):
    author_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:

        weight = relation['popularity'] * relation['possibility']
        risk = relation['popularity'] * relation['possibility'] * (1 - relation['possibility'])
        if risk == 0:
            continue
        if relation['actors'][0] not in author_dict:
            author_dict[relation['actors'][0]] = RiskNode(2)
        if relation['actors'][1] not in author_dict:
            author_dict[relation['actors'][1]] = RiskNode(2)
        author_dict[relation['actors'][0]].increase_neighbor(relation['actors'][1], 0, weight)
        author_dict[relation['actors'][1]].increase_neighbor(relation['actors'][0], 0, weight)
        author_dict[relation['actors'][0]].set_neighbor_risk(relation['actors'][1], -risk)
        author_dict[relation['actors'][1]].set_neighbor_risk(relation['actors'][0], -risk)
    return author_dict


def process_PPI_file(file_path):
    protein_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:

        weight = relation['weight'] * relation['possibility']
        risk = relation['weight'] * relation['possibility'] * (1 - relation['possibility'])
        if risk == 0:
            continue
        if relation['protein'][0] not in protein_dict:
            protein_dict[relation['protein'][0]] = RiskNode(2)
        if relation['protein'][1] not in protein_dict:
            protein_dict[relation['protein'][1]] = RiskNode(2)
        protein_dict[relation['protein'][0]].increase_neighbor(relation['protein'][1], 0, weight)
        protein_dict[relation['protein'][1]].increase_neighbor(relation['protein'][0], 0, weight)
        protein_dict[relation['protein'][0]].set_neighbor_risk(relation['protein'][1], -risk)
        protein_dict[relation['protein'][1]].set_neighbor_risk(relation['protein'][0], -risk)
    return protein_dict


def process_ReIDC_file(file_path):
    node_dict = {}
    text_file = open(file_path)
    line = text_file.readline()
    while line:
        possibility = line.split(" ")[2]
        if possibility != '0.000000':
            node1 = line.split(" ")[0]
            node2 = line.split(" ")[1]
            weight = float(possibility)
            risk = float(possibility) * (1-float(possibility))
            if node1 not in node_dict:
                node_dict[node1] = RiskNode(2)
            if node2 not in node_dict:
                node_dict[node2] = RiskNode(2)
            node_dict[node1].increase_neighbor(node2, 0, weight)
            node_dict[node2].increase_neighbor(node1, 0, weight)
            node_dict[node1].set_neighbor_risk(node2, -risk)
            node_dict[node2].set_neighbor_risk(node1, -risk)
        line = text_file.readline()
    return node_dict


def risk_averse_peel(node_dict, C_list, rho_list, B_list, precision=0.001):
    result = {'risk': dict(), 'weight': dict(), 'size': dict(), 'subgraph': dict()}
    for key in result:
        for C in C_list:
            result[key][C] = dict()
            for rho in rho_list:
                result[key][C][rho] = dict()

    pos_count = 0
    n = node_dict.__len__()
    print("initially the graph will have " + str(n) + " nodes")
    lambda2 = 1
    for node in node_dict.keys():
        for neighbor in node_dict[node].neighbor_dict.keys():
            #         if 0 in node_dict[node].neighbor_dict[neighbor]:
            pos_count += node_dict[node].neighbor_dict[neighbor][0]
    for B in B_list:

        for rho in rho_list:

            for C in C_list:
                print("!!!!!")
                print("Parameters set as: rho = {0}, B = {1}, C = {2}.".format(str(rho), str(B), str(C)))
                print("!!!!!")
                lambda1 = rho * lambda2
                low_bound = 0
                # To speed up the process, we usually set up bound to 20, which is bigger than most max q, instead of
                # the possible highest value below.
                up_bound = 20
                #  up_bound = (pos_count + lambda1 * n) / lambda2

                accelerate_flag = True
                while True:
                    # peeling with edge = pos - q * neg, and find if there exist a subgraph whose density > q * lambda2 - lambda1
                    # first build fib heap based on q
                    if accelerate_flag:
                        q = low_bound + (up_bound - low_bound) / 2
                    else:
                        q = (up_bound + low_bound) / 2
                    node_dict_q = {}
                    total_C_degree = 0
                    total_positive_degree = 0
                    fib_heap = FibonacciHeap()
                    for node in node_dict.keys():
                        node_dict_q[node] = SimpleNode()
                        for neighbor in node_dict[node].neighbor_dict.keys():
                            C_temp_degree_each = 0
                            positive_degree_each = 0
                            # here we already store disabled interactions as negative values
                            C_temp_degree_each += B * q * node_dict[node].neighbor_dict[neighbor][1]
                            C_temp_degree_each += C * node_dict[node].neighbor_dict[neighbor][0]
                            positive_degree_each += node_dict[node].neighbor_dict[neighbor][0]
                            node_dict_q[node].increase_neighbor(neighbor, C_temp_degree_each, positive_degree_each)
                            # to avoid influence from loop
                            if node == neighbor:
                                total_C_degree += C_temp_degree_each
                                total_positive_degree += positive_degree_each
                        node_dict_q[node].fib_node = fib_heap.insert(node_dict_q[node].C_degree, node)
                        total_C_degree += node_dict_q[node].C_degree
                        total_positive_degree += node_dict_q[node].positive_degree
                    total_positive_degree = total_positive_degree / 2
                    total_C_degree = total_C_degree / 2
                    exist_flag, C_avg, pos_avg, S_size, subgraph = peeling(node_dict_q, total_C_degree,
                                                        total_positive_degree, fib_heap, q, B, C, lambda1, lambda2)
                    print("current q={0}, find subgraph meet constraint={1}, average positive degree = {2}, the size "
                          "of the subgraph is {3}".format( str(q), str(exist_flag), str(pos_avg), str(S_size)))
                    if exist_flag:
                        if accelerate_flag:
                            accelerate_flag = False
                        if q - low_bound < precision and up_bound - q < precision:
                            weight = pos_avg
                            risk = (C*pos_avg - C_avg)/(B*q)

                            print("~~~~~~~~~~~~~")
                            print(
                                "rho, C, result q, corresponding (max density)subgraph's size, density, average risk:")
                            print(rho, C, q, S_size, weight, risk)
                            print("~~~~~~~~~~~~~")

                            result['size'][C][rho][B] = S_size
                            result['risk'][C][rho][B] = risk
                            result['weight'][C][rho][B] = weight
                            result['subgraph'][C][rho][B] = subgraph
                            break
                        else:
                            low_bound = q
                    else:
                        up_bound = q

    return result['weight'], result['risk'], result['size'], result['subgraph']


if __name__ == "__main__":
    uncertain_file = input('Enter the uncertain graph name(tmdb/dblp/ppi/reidc): ')
    file_path = input('Enter the file path(default if you input nothing): ')
    if len(file_path) == 0:
        file_path = '../datasets/tmdb/tmdb_2017.json'

    node_dict = {}
    # peeling preprocess for Tmdb
    if uncertain_file == 'tmdb':
        node_dict = process_tmdb_file(file_path)

    # peeling preprocess for DBLP
    elif uncertain_file == 'dblp':
        node_dict = process_dblp_file(file_path)

    # peeling preprocess for PPI datasets
    elif uncertain_file == 'ppi':
        node_dict = process_PPI_file(file_path)

    # peeling preprocess for PPI datasets
    elif uncertain_file == 'reidc':
        node_dict = process_ReIDC_file(file_path)

    else:
        assert Exception('uncertain dataset file type not expected!')

    B_list = [0.5, 1, 2, 5]
    C_list = [1]
    rho_list = [0.5, 1, 2, 5]
    weight, risk, size, subgraphs = risk_averse_peel(node_dict, C_list, rho_list, B_list)
    print('Average weight dictionary with indices as parameters C, rho and B.')
    print(weight)
    print('Average risk dictionary with indices as parameters C, rho and B.')
    print(risk)
    print('Subgraph size dictionary with indices as parameters C, rho and B.')
    print(size)
    print('Result subgraph dictionary with indices as parameters C, rho and B. (only subgraph with less than 30 nodes)')
    print(subgraphs)
