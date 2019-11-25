from lib.fibheap import FibonacciHeap
from lib.SimpleNode import SimpleNode
import json, sys, os


def peeling(node_dict, total_C_degree, total_positive_degree, fib_heap, q, B, C, lambda1, lambda2,
            max_subgraph_output_size=1000):
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
            if len(node_dict) < max_subgraph_output_size:
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


def process_uncertain_file(file_path, ignore_list=None):
    if ignore_list is None:
        ignore_list = []
    node_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:

        weight = relation['weight'] * relation['possibility']
        risk = relation['weight'] * relation['possibility'] * (1 - relation['possibility'])
        if risk == 0 or relation['nodes'][0] in ignore_list or relation['nodes'][1] in ignore_list:
            continue
        if relation['nodes'][0] not in node_dict:
            node_dict[relation['nodes'][0]] = RiskNode(2)
        if relation['nodes'][1] not in node_dict:
            node_dict[relation['nodes'][1]] = RiskNode(2)
        node_dict[relation['nodes'][0]].increase_neighbor(relation['nodes'][1], 0, weight)
        node_dict[relation['nodes'][1]].increase_neighbor(relation['nodes'][0], 0, weight)
        node_dict[relation['nodes'][0]].set_neighbor_risk(relation['nodes'][1], -risk)
        node_dict[relation['nodes'][1]].set_neighbor_risk(relation['nodes'][0], -risk)
    return node_dict


def process_signed_file(file_path):
    node_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:
        risk = -relation['positive']
        weight = -relation['negative']
        # weight = relation['positive']
        # risk = relation['negative']

        if relation['nodes'][0] not in node_dict:
            node_dict[relation['nodes'][0]] = RiskNode(2)
        if relation['nodes'][1] not in node_dict:
            node_dict[relation['nodes'][1]] = RiskNode(2)
        node_dict[relation['nodes'][0]].increase_neighbor(relation['nodes'][1], 0, weight)
        node_dict[relation['nodes'][1]].increase_neighbor(relation['nodes'][0], 0, weight)
        node_dict[relation['nodes'][0]].set_neighbor_risk(relation['nodes'][1], risk)
        node_dict[relation['nodes'][1]].set_neighbor_risk(relation['nodes'][0], risk)
    return node_dict


def build_fib_heap(node_dict, B, C, q):
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
    return node_dict_q, total_positive_degree, total_C_degree, fib_heap


def risk_averse_peel(node_dict, C_list, rho_list, B_list, precision=0.1):
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
                up_bound = 100
                #  up_bound = (pos_count + lambda1 * n) / lambda2

                accelerate_flag = True
                while True:
                    # peeling with edge = pos - q * neg, and find if there exist a subgraph whose density > q * lambda2 - lambda1
                    # first build fib heap based on q
                    if accelerate_flag:
                        q = low_bound + (up_bound - low_bound) / 2
                    else:
                        q = (up_bound + low_bound) / 2
                    node_dict_q, total_positive_degree, total_C_degree, fib_heap = build_fib_heap(node_dict, B, C, q)
                    exist_flag, C_avg, pos_avg, S_size, subgraph = peeling(node_dict_q, total_C_degree,
                                                                           total_positive_degree, fib_heap, q, B, C,
                                                                           lambda1, lambda2)
                    print("current q={0}, find subgraph meet constraint={1}, average positive degree = {2}, the size "
                          "of the subgraph is {3}".format(str(q), str(exist_flag), str(pos_avg), str(S_size)))
                    if exist_flag:
                        if accelerate_flag:
                            accelerate_flag = False
                        if q - low_bound < precision and up_bound - q < precision:
                            weight = pos_avg
                            risk = (C * pos_avg - C_avg) / (B * q)

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
    uncertain_file = input('Enter the uncertain graph name(uncertain/dblp/ppi/reidc): ')
    file_path = input('Enter the file path(default if you input nothing): ')
    if len(file_path) == 0:
        file_path = '../datasets/tmdb/tmdb_2017.json'

    node_dict = {}
    # peeling preprocess for Tmdb
    if uncertain_file == 'dblp':
        node_dict = process_dblp_file(file_path)

    # peeling preprocess for PPI datasets
    elif uncertain_file == 'uncertain':
        node_dict = process_uncertain_file(file_path)

    # DSD on a graph several times to get "almost densest subgraphs"
    elif uncertain_file == 'multiDSD':
        epsilon = input("Enter epsilon to define 'almost densest subgraph':")
        weight_list = {}
        risk_list = {}
        size_list = {}
        subgraph_list = {}
        almost = {}
        opt = 0
        for file in os.listdir(file_path):
            ignore_list = []
            count = 0
            if file == ".DS_Store":
                continue
            weight_list[file] = {}
            risk_list[file] = {}
            size_list[file] = {}
            subgraph_list[file] = {}
            for i in range(100):
                node_dict = process_uncertain_file(file_path + "/" + file, ignore_list)
                B = 1
                C = 1
                rho = 1
                if len(node_dict) <= 1:
                    break
                weight, risk, size, subgraphs = risk_averse_peel(node_dict, [C], [rho], [B])
                weight_list[file][i] = (weight[C][rho][B])
                if i == 0:
                    opt = weight[C][rho][B]
                else:
                    if weight[C][rho][B] >= opt * epsilon:
                        count += 1
                    else:
                        break
                risk_list[file][i] = risk[C][rho][B]
                subgraph_list[file][i] = subgraphs[C][rho][B]
                size_list[file][i] = size[C][rho][B]
                ignore_list = ignore_list + subgraphs[C][rho][B]
            almost[file] = count
        print("Count of almost densest subgraph")
        print(almost)
        print('Average weight dictionary with indices as file names and DSD iteration.')
        print(weight_list)
        print('Average risk dictionary with indices as file names and DSD iteration.')
        print(risk_list)
        print('Subgraph size dictionary with indices as file names and DSD iteration.')
        print(size_list)
        print(
            'Result subgraph dictionary with indices as file names and DSD iteration. (subgraphs with size out of '
            'bound will output as None)')
        print(subgraph_list)
        sys.exit()

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
    print('Result subgraph dictionary with indices as parameters C, rho and B. (subgraphs with size out of '
          'bound will output as None)')
    print(subgraphs)
