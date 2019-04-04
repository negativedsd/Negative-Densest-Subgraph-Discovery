from lib.fibheap import FibonacciHeap
from lib.SimpleNode import SimpleNode
import json


def peeling(node_dict, total_C_degree, total_q_degree, fib_heap, q, lambda1, lambda2):
    n = node_dict.__len__()
    avg_degree = (float)(total_C_degree) / n
    q_avg = total_q_degree / n
    # outputs we want
    max_C_avg = avg_degree
    S_size = n

    for i in range(n - 1):

        if i % 50000 == 0:
            print(i, max_C_avg, q_avg)

        # find min node from graph (remove from heap)
        node_to_remove = fib_heap.extract_min().value
        for neighbor in node_dict[node_to_remove].neighbor_dict.keys():

            # get dictionary that has all edges between two nodes
            C_degree_loss = node_dict[node_to_remove].neighbor_dict[neighbor][0]
            node_dict[neighbor].Cdegree -= C_degree_loss
            q_degree_loss = node_dict[node_to_remove].neighbor_dict[neighbor][1]
            node_dict[neighbor].qdegree -= q_degree_loss

            # here the key can be actually increased
            if neighbor != node_to_remove:
                fib_heap.decrease_key(node_dict[neighbor].fib_node, node_dict[neighbor].Cdegree)
                del node_dict[neighbor].neighbor_dict[node_to_remove]
            total_C_degree -= C_degree_loss
            total_q_degree -= q_degree_loss

        del node_dict[node_to_remove]
        avg_degree = (float)(total_C_degree) / (n - i - 1)
        q_avg = total_q_degree / (n - i - 1)
        S_size = n - i - 1
        if q_avg > q * lambda2 - lambda1:
            if len(node_dict)<100:
                print(list(node_dict))
            return True, avg_degree, q_avg, S_size

    return False, max_C_avg, q_avg, S_size


class RiskNode:
    degree = None
    total_degree = None
    neighbor_dict = None
    papercount = None

    def __init__(self, n):
        self.degree = [0] * n
        self.neighbor_dict = {}
        self.total_degree = 0
        self.papercount = 0

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
        actor_dict[relation['actors'][0]].increase_neighbor(relation['actors'][1],0,weight)
        actor_dict[relation['actors'][1]].increase_neighbor(relation['actors'][0],0,weight)
        actor_dict[relation['actors'][0]].set_neighbor_risk(relation['actors'][1],-risk)
        actor_dict[relation['actors'][1]].set_neighbor_risk(relation['actors'][0],-risk)
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
        author_dict[relation['actors'][0]].increase_neighbor(relation['actors'][1],0,weight)
        author_dict[relation['actors'][1]].increase_neighbor(relation['actors'][0],0,weight)
        author_dict[relation['actors'][0]].set_neighbor_risk(relation['actors'][1],-risk)
        author_dict[relation['actors'][1]].set_neighbor_risk(relation['actors'][0],-risk)
    return author_dict

def process_PPI_file(file_path):
    author_dict = {}
    relation_list = json.load(open(file_path))
    for relation in relation_list:

        weight = relation['weight'] * relation['possibility']
        risk = relation['weight'] * relation['possibility'] * (1 - relation['possibility'])
        if risk == 0:
            continue
        if relation['protein'][0] not in author_dict:
            author_dict[relation['protein'][0]] = RiskNode(2)
        if relation['protein'][1] not in author_dict:
            author_dict[relation['protein'][1]] = RiskNode(2)
        author_dict[relation['protein'][0]].increase_neighbor(relation['protein'][1],0,weight)
        author_dict[relation['protein'][1]].increase_neighbor(relation['protein'][0],0,weight)
        author_dict[relation['protein'][0]].set_neighbor_risk(relation['protein'][1],-risk)
        author_dict[relation['protein'][1]].set_neighbor_risk(relation['protein'][0],-risk)
    return author_dict


def risk_averse_peel(uncertain_file='tmdb', filepath='../datasets/tmdb/tmdb_2017.json', precision=0.02):
    # peeling preprocess for Tmdb
    if uncertain_file == 'tmdb':
        node_dict = process_tmdb_file(filepath)

    # peeling preprocess for DBLP
    elif uncertain_file == 'dblp':
        node_dict = process_dblp_file(filepath)

    # peeling preprocess for PPI datasets
    elif uncertain_file == 'ppi':
        node_dict = process_PPI_file(filepath)

    else:
        assert Exception('uncertain dataset file type not expected!')

    result = {'risk':{},'weight':{},'size':{}}

    pos_count = 0
    n = node_dict.__len__()
    print("initially the graph will have " + str(n) + " nodes")
    lambda2 = 1
    for node in node_dict.keys():
        for neighbor in node_dict[node].neighbor_dict.keys():
            #         if 0 in node_dict[node].neighbor_dict[neighbor]:
            pos_count += node_dict[node].neighbor_dict[neighbor][0]
    rho_list = [0.5, 1, 2, 10]
    C_list = [0.25,0.5,1,2,3,4,5,6]
    for rho in rho_list:
        print("!!!!!")
        print("now we have rho as ", rho)
        print("!!!!!!")
        for C in C_list:
            print("!!!!!")
            print("now we have C as ", C)
            print("!!!!!!")
            lambda1 = rho * lambda2
            lowbound = 0
            upbound = 20
            #  upbound = (pos_count + lambda1 * n) / lambda2

            accerate_flag = True
            while True:
                #     peeling with edge = pos - q * neg, and find if there exist a subgraph whose density > q * lambda2 - lambda1
                # first build fib heap based on q
                if accerate_flag:
                    q = lowbound + (upbound - lowbound) / 2
                else:
                    q = (upbound + lowbound) / 2
                node_dict_q = {}
                total_C_degree = 0
                total_q_degree = 0
                fib_heap = FibonacciHeap()
                for node in node_dict.keys():
                    node_dict_q[node] = SimpleNode()
                    for neighbor in node_dict[node].neighbor_dict.keys():
                        C_temp_degree_each = 0
                        q_temp_degree_each = 0
                        # here we already store disabled interactions as negative values
                        C_temp_degree_each += node_dict[node].neighbor_dict[neighbor][1]
                        q_temp_degree_each += q * node_dict[node].neighbor_dict[neighbor][1]
                        C_temp_degree_each += C * node_dict[node].neighbor_dict[neighbor][0]
                        q_temp_degree_each += node_dict[node].neighbor_dict[neighbor][0]
                        node_dict_q[node].increase_neighbor(neighbor, C_temp_degree_each, q_temp_degree_each)
                        # to avoid influence from loop
                        if node == neighbor:
                            total_C_degree += C_temp_degree_each
                            total_q_degree += q_temp_degree_each
                    node_dict_q[node].fib_node = fib_heap.insert(node_dict_q[node].Cdegree, node)
                    total_C_degree += node_dict_q[node].Cdegree
                    total_q_degree += node_dict_q[node].qdegree
                #  print(total_C_degree, total_q_degree)
                total_q_degree = total_q_degree / 2
                total_C_degree = total_C_degree / 2
                print(total_C_degree, total_q_degree)
                exist_flag, max_avg, q_avg, S_size = peeling(node_dict_q, total_C_degree, total_q_degree, fib_heap, q,
                                                             lambda1, lambda2)
                print(q, exist_flag, max_avg, q_avg, S_size)
                if exist_flag:
                    if accerate_flag:
                        accerate_flag = False
                    if q - lowbound < precision and upbound - q < precision:
                        weight = (max_avg-q_avg/q)/(C-1/q)
                        risk = C*weight - max_avg
                        print("~~~~~~~~~~~~~")
                        print("rho, C, result q, corresponding (max density)subgraph's size, density, average risk:")
                        print(rho, C, q, S_size, weight, risk)
                        if C not in result['size']:
                            result['size'][C] = {}
                        result['size'][C][rho] = S_size
                        if C not in result['risk']:
                            result['risk'][C] = {}
                        result['risk'][C][rho] = risk
                        if C not in result['weight']:
                            result['weight'][C] = {}
                        result['weight'][C][rho] = weight
                        print("~~~~~~~~~~~~~")
                        break
                    else:
                        lowbound = q
                else:
                    upbound = q

    print(result)


if __name__ == "__main__":
    risk_averse_peel()
