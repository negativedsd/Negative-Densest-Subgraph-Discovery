from lib.fibheap import FibonacciHeap
import copy


class Node:
    degree = None
    # {'name of neighbor':{'edge type':weight}}
    neighbor_dict = None
    fib_node = None

    def __init__(self):
        # self.name = name
        self.degree = 0
        self.neighbor_dict = {}

    def increase_neighbor(self, name, type, degree):
        if name not in self.neighbor_dict:
            self.neighbor_dict[name] = {type: degree}
        else:
            if type not in self.neighbor_dict[name]:
                self.neighbor_dict[name][type] = degree
            else:
                self.neighbor_dict[name][type] += degree
        self.degree += degree


def read_graph_file(path_list, t_list, enabled_list, interaction_degree, neg_value, C, pos_value_one):
    node_dict = {}
    fib_heap = FibonacciHeap()
    total_degree = 0
    edge_count_dict = {}
    for i in range(len(path_list)):
        print('reading file: ', path_list[i])
        with open(path_list[i]) as file:
            line = file.readline()
            count = 0
            while line:
                if 'deleted' in line:
                    line = file.readline()
                    continue
                count += 1
                if count % 100000 == 0:
                    print(count)
                line_list = [(s) for s in line.split() if s.isdigit()]
                try:
                    if pos_value_one:
                        pos = 1
                    else:
                        pos = (int)(line_list[2])
                except ValueError:
                    print(line)
                    line = file.readline()
                    continue
                if enabled_list[i] == 0:
                    pos = 0
                line = file.readline()
                if enabled_list[i] == 1:
                    neg = 0
                else:
                    neg = neg_value
                edge_type = t_list[i]
                if edge_type not in edge_count_dict:
                    edge_count_dict[edge_type] = 0
                edge_count_dict[edge_type] += 1
                edge_degree = C * pos - neg
                if edge_type not in interaction_degree:
                    interaction_degree[edge_type] = 0
                interaction_degree[edge_type] += edge_degree
                if line_list[0] not in node_dict:
                    node_dict[line_list[0]] = Node()
                node_dict[line_list[0]].increase_neighbor(line_list[1], edge_type, edge_degree)
                if line_list[0] != line_list[1]:
                    if line_list[1] not in node_dict:
                        node_dict[line_list[1]] = Node()

                    node_dict[line_list[1]].increase_neighbor(line_list[0], edge_type, edge_degree)

                total_degree += edge_degree

            print('read complete', count, 'lines')

    print('building Fib heap')
    # check_degree = 0
    # for key in node_dict.keys():
    #     node_dict[key].fib_node = fib_heap.insert(node_dict[key].degree, key)
    #     check_degree += node_dict[key].degree
    # print(total_degree, check_degree)
    return fib_heap, node_dict, total_degree, edge_count_dict


def get_sum(dictionary):
    result = 0
    for key in dictionary.keys():
        result += dictionary[key]
    return result


def get_min_index(l):
    minIndex = 0
    minValue = l[0][0]
    for i in range(1, len(l)):
        if minValue > l[i][0]:
            minIndex = i
    return minIndex


def peeling(node_dict, total_degree, interaction_degree, edge_count_dict, fib_heap):
    n = node_dict.__len__()
    avg_degree = total_degree / n
    result = []
    # outputs we want
    max_avg = avg_degree
    S_size = n
    max_interaction_degree = copy.deepcopy(interaction_degree)
    for key in max_interaction_degree:
        max_interaction_degree[key] = interaction_degree[key] / S_size

    result.append([max_avg, S_size, max_interaction_degree])
    min_index = 0

    print(n, 'nodes')
    for i in range(n - 1):
        if i % 100000 == 0:
            print(i)
            print(result)
        # find min node from graph (remove from heap)
        node_to_remove = fib_heap.extract_min().value

        # for every neighbor node this min node have
        for neighbor in node_dict[node_to_remove].neighbor_dict.keys():

            # get dictionary that has all edges between two nodes
            one_neighbor_edge_dict = node_dict[node_to_remove].neighbor_dict[neighbor]

            degree_loss = 0
            for edge_type in one_neighbor_edge_dict:
                if edge_type in edge_count_dict:
                    interaction_degree[edge_type] -= one_neighbor_edge_dict[edge_type]
                    # edge_count_dict[edge_type] -= 1
                else:
                    print('!!!')
                degree_loss += one_neighbor_edge_dict[edge_type]

            node_dict[neighbor].degree -= degree_loss

            # here the key can be actually increased
            if neighbor != node_to_remove:
                fib_heap.decrease_key(node_dict[neighbor].fib_node, node_dict[neighbor].degree)
                del node_dict[neighbor].neighbor_dict[node_to_remove]
            total_degree -= degree_loss

        del node_dict[node_to_remove]
        avg_degree = total_degree / (n - i - 1)
        #     print(avg_degree)
        if result[min_index][0] < avg_degree:
            max_avg = avg_degree
            S_size = n - i - 1
            max_interaction_degree = copy.deepcopy(interaction_degree)
            for key in max_interaction_degree:
                max_interaction_degree[key] = interaction_degree[key] / S_size

            if len(result) < 5:
                result.append([max_avg, S_size, max_interaction_degree])
            else:
                result[min_index] = [max_avg, S_size, max_interaction_degree]
                min_index = get_min_index(result)
    return result


def get_densest_subgraph(interactions=None, neg_value=1, C=1):
    if interactions != None:
        interaction_list = interactions
    else:
        interaction_list = ['retweet', 'reply']
    interaction_degree = {}
    results = {0: {}, 1: {}, 2: {}}
    enable_list = [[0, 1], [1, 0], [1, 1]]
    for day in range(1, 8):
        for index in range(len(enable_list)):
            enable_vector = enable_list[index]
            for i in range(len(enable_vector)):
                interaction_degree[interaction_list[i]] = 0
            date_form = '-2018-02-'
            if day < 10:
                date_form += '0'
            date_form += str(day)
            fib_heap, node_dict, total_degree, edge_count_dict = read_graph_file(
                [interaction_list[0] + date_form + '.txt',
                 interaction_list[1] + date_form + '.txt'],
                interaction_list, enable_vector,
                interaction_degree, neg_value,
                C, True)
            results[index][day] = peeling(node_dict, total_degree, interaction_degree, edge_count_dict, fib_heap)
            print('peeling for day', day, 'with enable list', enable_list[index], 'finished')
    print('!!!!!!!!!')
    print("The results (top 5 densest subgraph everyday in format [density, size, density for each interaction]):")
    print(results)
    print('!!!!!!!!!')

    max_degree_dict = [{}, {}, {}]
    for i in range(0, 3):
        for day in range(1, 8):
            temp_max = 0
            for each in results[i][day]:
                if each[0] > temp_max:
                    max_degree_dict[i][day] = each
                    temp_max = each[0]
    print("The results (densest subgraph everyday in format [density, size, density for each interaction]):")
    print(max_degree_dict)


if __name__ == "__main__":
    get_densest_subgraph()
