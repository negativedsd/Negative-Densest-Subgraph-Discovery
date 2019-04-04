from lib.fibheap import FibonacciHeap
from lib.SimpleNode import SimpleNode
from synthetic.erdos import *
import pickle
import matplotlib.pyplot as plt
import numpy as np


def peeling(node_dict, total_C_degree, total_q_degree, fib_heap, q, lambda1, lambda2):
    n = node_dict.__len__()
    avg_degree = (float)(total_C_degree) / n
    q_avg = total_q_degree / n
    # outputs we want
    max_C_avg = avg_degree
    S_size = n

    for i in range(n - 1):

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
        # if max_C_avg < avg_degree:
        max_C_avg = avg_degree
        q_avg = total_q_degree / (n - i - 1)
        S_size = n - i - 1
        if q_avg > q * lambda2 - lambda1:

            return True, max_C_avg, q_avg, S_size, node_dict.keys()
    return False, max_C_avg, q_avg, S_size, node_dict.keys()


def binarySearch(pos_count, lambda1, lambda2, C, node_dict, precision = 0.01):
    lowbound = 0
    n = node_dict.__len__()
    upbound = (pos_count + lambda1 * n) / lambda2
    print('initial upbound:',upbound)
    while True:
        #     peeling with edge = C * pos - neg, and find if there exist a subgraph whose density > q * lambda2 - lambda1
        # first build fib heap based on q
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
                for edge in node_dict[node].neighbor_dict[neighbor]:
                    if node_dict[node].neighbor_dict[neighbor][edge] < 0:
                        #                         here we already store disabled interactions as negative values
                        C_temp_degree_each += node_dict[node].neighbor_dict[neighbor][edge]
                        q_temp_degree_each += q * node_dict[node].neighbor_dict[neighbor][edge]
                    else:
                        C_temp_degree_each += C * node_dict[node].neighbor_dict[neighbor][edge]
                        q_temp_degree_each += node_dict[node].neighbor_dict[neighbor][edge]
                node_dict_q[node].increase_neighbor(neighbor, C_temp_degree_each, q_temp_degree_each)
                # to avoid influence from self connect edges
                if node == neighbor:
                    total_C_degree += C_temp_degree_each
                    total_q_degree += q_temp_degree_each
            node_dict_q[node].fib_node = fib_heap.insert(node_dict_q[node].Cdegree, node)
            total_C_degree += node_dict_q[node].Cdegree
            total_q_degree += node_dict_q[node].qdegree
        total_q_degree = total_q_degree / 2
        total_C_degree = total_C_degree / 2
        # print(total_C_degree, total_q_degree)
        exist_flag, max_avg, q_avg, S_size, output_nodes = peeling(node_dict_q, total_C_degree, total_q_degree, fib_heap, q, lambda1, lambda2)
        # print(q, exist_flag, max_avg, q_avg, S_size, q * lambda2 - lambda1)
        if exist_flag:
            if q - lowbound < precision:
                # print("~~~~~~~~~~~~~")
                print("result q, corresponding (max density)subgraph's size:", q, S_size)
                print('the nodes are:', output_nodes)
#                 temp_list.append(S_size)
                print("~~~~~~~~~~~~~")
                return S_size
            else:
                lowbound = q
        else:
            upbound = q


def countPositiveEdge(node_dict):
    pos_count = 0
    for node in node_dict.keys():
        for neighbor in node_dict[node].neighbor_dict.keys():
            if 0 in node_dict[node].neighbor_dict[neighbor]:
                pos_count += 1
    return pos_count


def queryWithDiffRhoAndC(rho_list, C_list, graph_file_path=None, ratio=1.0, lambda2=1.0):

#     :param USE_GOOD_GRAPH: True to use a pickled graph that shows with different C the algorithm can get different planted part,
#     False to generate a new random graph
#     :param ratio: Expected weight of positive edges / Expected weight of negative edges
#     :param rho_list: list of rho value used to construct lambda pairs
#     :param C_list: list of C value used in greedy peeling algorithm
#     :param lambda2: value of lambda2, and lambda1 = rho * lambda2
#     :return: a list, each element corresponds to a different C value, and the element itself is also a
#     list, each element corresponds to a different rho value, and the element itself is the output sizes

    result = []

    if graph_file_path == None:
        node_dict = erdos_renyi(200, 0.02, ratio)
        node_dict = plant(node_dict, 12, 10, 20, 5, False)
    else:
        node_dict = pickle.load( open( graph_file_path, "rb" ) )

    n = node_dict.__len__()
    print("initially the graph will have " + str(n) + " nodes")

    pos_count = countPositiveEdge(node_dict)

    for C in C_list:
        print("!!!!!")
        print("now we have C as ", C)
        # print("!!!!!!")
        temp_list = []
        for rho in rho_list:
            lambda1 = rho * lambda2
            print("!!!!!")
            print("now we have lambda1 as ", lambda1, "lambda2 as ", lambda2)
            # print("!!!!!!")
            temp_list.append(binarySearch(pos_count, lambda1, lambda2, C, node_dict))
        result.append(temp_list)
    return result,node_dict


def queryWithDiffLambdas(lambda1_list, lambda2_list, graph_file_path=None, ratio=1.0, C=1.0):
    result = []
    if graph_file_path == None:
        node_dict = erdos_renyi(200, 0.01, ratio)
        node_dict = plant(node_dict, 12, 10, 20, 6, False)
    else:
        node_dict = pickle.load( open( graph_file_path, "rb" ) )

    n = node_dict.__len__()
    print("initially the graph will have " + str(n) + " nodes")

    pos_count = countPositiveEdge(node_dict)

    for lambda1 in lambda1_list:
        temp_list = []
        for lambda2 in lambda2_list:
            print("!!!!!")
            print("now we have lambda1 as ", lambda1, "lambda2 as ", lambda2)
            # print("!!!!!!")
            temp_list.append(binarySearch(pos_count, lambda1, lambda2, C, node_dict))
        result.append(temp_list)
    return result,node_dict


def queryWithDiffPlantSize(size1_list, size2_list, lambda1=1, lambda2=1, ratio=1.0, C=1.0):
    result = []
    for size1 in size1_list:
        temp_list = []
        for size2 in size2_list:
            if size2<size1:
                continue
            node_dict = erdos_renyi(200, 0.01, ratio)
            node_dict = plant(node_dict, size1, 6, size2, 6, False)
            n = node_dict.__len__()

            sub1_degree = 0
            sub1_pos_degree = 0
            for i in range(0,size1):
                for j in range(0,size1):
                    if j in node_dict[i].neighbor_dict:
                        if 0 in node_dict[i].neighbor_dict[j]:
                            sub1_degree += node_dict[i].neighbor_dict[j][0]
                            sub1_pos_degree += node_dict[i].neighbor_dict[j][0]
                        if 1 in node_dict[i].neighbor_dict[j]:
                            sub1_degree += node_dict[i].neighbor_dict[j][1]
            print('density of subgraph1 is:', sub1_degree/size1)
            print('Positive density of subgraph1 is:', sub1_pos_degree / size2)

            sub2_degree = 0
            sub2_pos_degree = 0
            for i in range(200-size2, 200):
                for j in range(200-size2, 200):
                    if j in node_dict[i].neighbor_dict:
                        if 0 in node_dict[i].neighbor_dict[j]:
                            sub2_degree += node_dict[i].neighbor_dict[j][0]
                            sub2_pos_degree += node_dict[i].neighbor_dict[j][0]
                        if 1 in node_dict[i].neighbor_dict[j]:
                            sub2_degree += node_dict[i].neighbor_dict[j][1]
            print('density of subgraph2 is:', sub2_degree/size2)
            print('Positive density of subgraph2 is:', sub2_pos_degree / size2)

            print("initially the graph will have " + str(n) + " nodes")

            pos_count = countPositiveEdge(node_dict)

            print("!!!!!")
            print("now we have size1 as ", size1, "size2 as ", size2)
            temp_list.append(binarySearch(pos_count, lambda1, lambda2, C, node_dict))
        result.append(temp_list)
    return result


# print('Experiment with different rho and C')
# rho_list = [0, 1, 5, 10,30,1000]
# C_list = [0.1,1,10,100]
# result1, node_dict1 = queryWithDiffRhoAndC(rho_list, C_list, 'good_rho_graph.p')
# print('Output size list used for plotting:', result1)
# print('\n\n')

# print('Experiment with different lambda1 and lambda2')
# lambda1_list = [0.1, 1, 10,50,100]
# lambda2_list = [0.1, 1, 10,50,100]
# result2, node_dict2 = queryWithDiffLambdas(lambda1_list,lambda2_list, graph_file_path='good_C_graph.p', C=1)
# print('Output size list used for plotting:', result2)

print('Experiment with different lambda1 and lambda2')
size1_list = [9, 12, 15, 20]
size2_list = [9, 12, 15, 20]
result3 = queryWithDiffPlantSize(size1_list, size2_list, lambda1=10, lambda2=1)
print('Output size list used for plotting:', result3)
