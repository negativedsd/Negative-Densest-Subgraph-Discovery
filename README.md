# Novel Dense Subgraph Discovery Primitives: Risk Aversion and Exclusion Queries 

## Motivation

In the densest subgraph problem, given a weighted undirected graph G(V,E,w) with non-negative edge weights, we are asked to find a subset of nodes S that maximizes the degree density w(S)/|S|, where w(S) is the sum of the edge weights induced by SS. This problem is a well studied problem, known as the  *densest subgraph problem*, and is solvable in polynomial time. But what happens when the edge weights are *negative*? Is the problem still solvable in polynomial time? Also, why should we care about the densest subgraph problem in the presence of negative weights?  

In this work we answer the aforementioned question. Specifically, we provide two novel graph mining primitives that are applicable to a wide variety of applications. Our primitives can be used to answer questions such as ``how can we find a dense subgraph in Twitter with lots of replies and mentions but no follows?'', ``how do we extract a dense subgraph with high expected reward and low risk from an uncertain graph''? We formulate both problems mathematically as special instances of dense subgraph discovery in graphs with negative weights.  

## Code 

For a description on how to use our code, and a demo, check the Jupyter notebook demo.ipynb.

## Authors 

- [Charalampos E. Tsourakakis](https://tsourakakis.com/)
- Tianyi Chen 
- Naonori Nakimura 
- Jakub Pachocki 
