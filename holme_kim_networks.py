#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:11:35 2022

@author: pkundu
"""

import networkx as nx

nn=300; mm=5; p=0.5; 
#gs=nx

for i in range(1,101):
    G=nx.powerlaw_cluster_graph(nn,mm,p)
#edgelist=nx.write_edgelist(G,"test.edgelist")
    nx.write_edgelist(G, 'E' + str(i) + '.txt', delimiter=' ', data = False)
#nx.savetxt(edgelist)
#nx.plot(G)
# nx.savetxt()
    