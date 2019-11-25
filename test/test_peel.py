import unittest
from lib.fibheap import FibonacciHeap
from risk_averse.peel_risk_averse import *


def build_graph():
	node_a = RiskNode(2)
	node_b = RiskNode(2)
	node_c = RiskNode(2)
	node_d = RiskNode(2)
	node_a.increase_neighbor('b', 0, 0.2)
	node_a.set_neighbor_risk('b', -(2 * 0.1 * (1 - 0.1)))
	node_a.increase_neighbor('c', 0, 0.2)
	node_a.set_neighbor_risk('c', -(1 * 0.2 * (1 - 0.2)))
	node_b.increase_neighbor('a', 0, 0.2)
	node_b.set_neighbor_risk('a', -(2 * 0.1 * (1 - 0.1)))
	node_c.increase_neighbor('a', 0, 0.2)
	node_c.set_neighbor_risk('a', -(1 * 0.2 * (1 - 0.2)))

	node_b.increase_neighbor('c', 0, 0.2)
	node_b.set_neighbor_risk('c', -(1 * 0.2 * (1 - 0.2)))
	node_c.increase_neighbor('b', 0, 0.2)
	node_c.set_neighbor_risk('b', -(1 * 0.2 * (1 - 0.2)))

	node_b.increase_neighbor('d', 0, 0.6)
	node_b.set_neighbor_risk('d', -(2 * 0.3 * (1 - 0.3)))
	node_d.increase_neighbor('b', 0, 0.6)
	node_d.set_neighbor_risk('b', -(2 * 0.3 * (1 - 0.3)))
	return {'a': node_a, 'b': node_b, 'c': node_c, 'd': node_d}


def isclose(a, b):
	return abs(a - b) <= 1e-05


class MyTestCase(unittest.TestCase):

	def test_json_load(self):
		node_dict = process_uncertain_file("./test.json")
		build_dict = build_graph()
		self.assertDictEqual(node_dict['a'].neighbor_dict, build_dict['a'].neighbor_dict)
		self.assertDictEqual(node_dict['b'].neighbor_dict, build_dict['b'].neighbor_dict)
		self.assertDictEqual(node_dict['c'].neighbor_dict, build_dict['c'].neighbor_dict)
		self.assertDictEqual(node_dict['d'].neighbor_dict, build_dict['d'].neighbor_dict)

	def test_single_peel_empty_result(self):
		build_dict = process_uncertain_file("./test.json")
		B = 1
		C = 1
		q = 1
		lambda1 = 1
		lambda2 = 10
		node_dict_q, total_positive_degree, total_C_degree, fib_heap = build_fib_heap(build_dict, B, C, q)
		exist_flag, C_avg, pos_avg, S_size, subgraph = peeling(node_dict_q, total_C_degree,
															   total_positive_degree, fib_heap, q, B, C, lambda1,
															   lambda2)

		self.assertEqual(exist_flag, False)
		self.assertEqual(isclose(C_avg, 0), True)
		self.assertEqual(isclose(pos_avg, 0), True)
		self.assertEqual(S_size, 1)

	def test_single_peel_normal(self):
		build_dict = process_uncertain_file("./test.json")
		B = 1
		C = 1
		q = 1
		lambda1 = 1
		lambda2 = 1.07
		node_dict_q, total_positive_degree, total_C_degree, fib_heap = build_fib_heap(build_dict, B, C, q)
		exist_flag, C_avg, pos_avg, S_size, subgraph = peeling(node_dict_q, total_C_degree,
															   total_positive_degree, fib_heap, q, B, C, lambda1,
															   lambda2)
		print(exist_flag, C_avg, pos_avg, S_size, subgraph)
		self.assertEqual(exist_flag, True)
		self.assertEqual(isclose(C_avg, 0.22 / 3), True)
		self.assertEqual(isclose(pos_avg, 0.8 / 3), True)
		self.assertEqual(S_size, 3)

	def test_risk_averse(self):
		build_dict = process_uncertain_file("./test.json")
		B_list = [1]
		C_list = [1]
		rho_list = [1]
		weight, risk, size, subgraphs = risk_averse_peel(build_dict, C_list, rho_list, B_list)
		truth_result_weight = 0.3
		truth_result_risk = 0.21
		self.assertEqual(isclose(weight[1][1][1], truth_result_weight), True)
		self.assertEqual(isclose(risk[1][1][1], truth_result_risk), True)
		self.assertEqual(size[1][1][1], 2)
		self.assertEqual(subgraphs[1][1][1], ['b','d'])


if __name__ == '__main__':
	unittest.main()
