import networkx as nx
import matplotlib.pyplot as plt 
import random 
import pandas as pd
from sys import argv

# The partial duplication model 
def partial_dupli(s, q, p, N):
	G=nx.Graph()
	G.add_nodes_from(range(s))
	for v in range(s, N):
		G.add_node(v)
		u = random.randint(0,max(G.nodes())-1)
		ns = list(nx.all_neighbors(G, u))
		if random.random() < q: #uniformally
			G.add_edge(u,v)
		for n in ns:
			if random.random() < p:
				G.add_edge(v,n)
	return G

# The duplication-divergence model
def duplication_div(s, q, p, r, N):
	ns = []
	G=nx.Graph()
	G.add_nodes_from(range(s))
	for v in range(s, N):
		G.add_node(v)
		u = random.randint(0,max(G.nodes())-1)
		ns = list(nx.all_neighbors(G, u))
		if random.random() < q:
			G.add_edge(u,v)
		for n in ns:
			G.add_edge(v,n)
			if random.random() > p: 
				if r < random.random():
					G.remove_edge(u,n)
				else:
					G.remove_edge(v,n)	
	return G

#Confirm that the produces networks which follow the Power law by plotting the degree distribution on a log-log scale.
#Producing the divergence duplication model:
G_dupli = duplication_div(1, 0.2, 0.3, 0.5, 1000)
plt.figure()
nx.draw(G_dupli,  node_size = 10)
plt.savefig("Graph_dupli.png", format="PNG", node_color='blue')

#Producing the divergence duplication model:
G_part = partial_dupli(1, 0.2, 0.3, 1000)
plt.figure()
nx.draw(G_part,  node_size = 10)
plt.savefig("Graph_part.png", format="PNG", node_color='blue')

#Plotting the Degree Distribution of the 2 graphs
degree_sequence= nx.degree_histogram(G_dupli) 
plt.figure()
plt.loglog(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot of duplication divergence model")
plt.ylabel("degree")
plt.xlabel("rank")
plt.savefig('degree_rank_duplication')

degree_sequence= nx.degree_histogram(G_part)
plt.figure()
plt.loglog(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot of partial model")
plt.ylabel("degree")
plt.xlabel("rank")
plt.savefig('degree_rank_partial_model')

# Compute the clustering coefficient for each network. Is the clustering coefficient maintained as the networks become larger?

print 'The average clustering coefficient of the model divergence duplication is:', nx.average_clustering(G_dupli)
print 'The average clustering coefficient of the model partial duplication is:', nx.average_clustering(G_part)

# Is the clustering coefficient maintained as the networks become larger?
for i in range(3000, 7000, 1000):
	G_dupli = partial_dupli(1, 0.2, 0.3, i)
	print 'The number of nodes of the network is: %d' % (i)
	print 'The average_clustering coefficient is:', nx.average_clustering(G_dupli)

# For the Yeast protein interaction network:
yeast = pd.read_csv(argv[1], sep='\s')
print ' The head of dataframe before filtering:\n'
print yeast.head()
print ' The head of dataframe after filtering:\n'
yeast_sub = yeast.loc[(yeast.combined_score >= 900)] 

print yeast_sub.head()

points = zip(yeast_sub.protein1, yeast_sub.protein2)
G_yeast=nx.Graph()
G_yeast.add_edges_from(points)
degree_sequence_yeast= nx.degree_histogram(G_yeast) # degree sequence
plt.figure()
plt.loglog(degree_sequence_yeast,'b-',marker='o')
plt.title("Degree rank plot of yeast data")
plt.ylabel("degree")
plt.xlabel("rank")

# Making a function to calculate average shortest pathway and giving a restriction for unconnected graphs
def evaluator(G):
	calc = list()
	ev1 = nx.average_clustering(G)
	if nx.is_connected(G) == True:
		ev2 = nx.average_shortest_path_length(G)
	else:
		for sub in nx.connected_component_subgraphs(G):
			if len(sub.nodes()) > 1:
				calc.append(nx.average_shortest_path_length(sub))
		ev2 = sum(calc)/len(calc)
	print 'Average clustering and average shortest path length coefficients:', (ev1, ev2)

# Compute the clustering coefficient.
# Compute the average shortest path-length.
print 'The topological profile of your data is:'
evaluator(G_yeast)
# These are the average clustering and average shortest path length coefficients of the Yeast graph before filtering (took 1 hour)
#(0.2775705795995122, 1.5528623151785386)

# Now generate an Erdos-Renyi graph, a partial duplication network and a duplication-divergence network with the same number of nodes as the yeast PIN.

nodes = nx.number_of_nodes(G_yeast)
print 'The number of nodes in the given network is:', nodes
print 'The number of edges in the the given network is:', G_yeast.number_of_edges()

# Plot the degree distribution an of each network, including the yeast PIN. Do they look as expected? Describe how and why the degree distributions look different.

G_random = nx.fast_gnp_random_graph(nodes, 0.01)
degree_sequence_random = nx.degree_histogram(G_random)

G_partial_dupli = partial_dupli(1, 0.2, 0.3, nodes)
degree_sequence_partial = nx.degree_histogram(G_partial_dupli) 

G_dupli_div= duplication_div(1, 0.2, 0.3, 0.5, nodes)
degree_sequence_div = nx.degree_histogram(G_dupli_div) 

fig = plt.figure()
ax1 = fig.add_subplot(411)
ax1.loglog(degree_sequence_partial,'b-',marker='o')
plt.title("Degree rank of partial duplication network")
ax2 = fig.add_subplot(412)
ax2.loglog(degree_sequence_div,'b-',marker='o')
plt.title("Degree rank of duplication divergence network")
ax3 = fig.add_subplot(413)
ax3.loglog(degree_sequence_random,'b-',marker='o')
plt.title("Degree rank of random network")
ax4 = fig.add_subplot(414)
ax4.loglog(degree_sequence_yeast,'b-',marker='o')
plt.title("Degree rank of Yeast network")
plt.savefig('comparison_degree.png')


# For each of the four networks (Erdos-Renyi, partial duplication, duplication-divergence and yeast PIN), try to remove random vertices. That is, for a given network, pick a random vertex and remove it from the network. Repeat this process and compute the degree distribution at regular intervals (for example, after removing 100, 200, ... vertices). This simulates a random attack on the network. Explain what happens to the degree distribution as the network is attacked. What happens to the clustering coefficient? The average shortest path-length?
#########################
# Making random attacks
#########################
# Random attack:
def random_attack(G, n):
	nodes_list = G.nodes()
	my_randoms = random.sample(nodes_list, n) #Without replacement
	G.remove_nodes_from(my_randoms)
	return G

# random network
print '\nRandom attack in randomized network\n'
G_random = nx.fast_gnp_random_graph(nodes, 0.01)
for i in range(1, 4): 
	G_sub = random_attack(G_random, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_randomized_random' + str(i*100) +'.png', format ='png')  
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o') #log(y/x)
		plt.title("Degree rank plot of random model - random attack")
		plt.savefig('random_network_%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# yeast network
print '\nRandom attack in yeast_network \n'
for i in range(1, 4): 
	G_sub = random_attack(G_yeast, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_yeast_random' + str(i*100) +'.png', format ='png')  
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of yeast model - random attack")
		plt.savefig('yeast_network_%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# partial duplication
print '\nRandom attack in partial duplication network \n'
G_partial_dupli = partial_dupli(1, 0.2, 0.3, nodes)
for i in range(1, 4): 
	G_sub = random_attack(G_partial_dupli, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_part_dupli_random' + str(i*100) +'.png', format ='png')  
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of partial duplication model - random attack")
		plt.savefig('partial_duplication_network_%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# divergence duplication
print '\nRandom attack in divergence duplication network \n'
G_dupli_div = duplication_div(1, 0.2, 0.3, 0.5, nodes)
for i in range(1, 4): 
	G_sub = random_attack(G_dupli_div, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_div_dupli_random' + str(i*100) +'.png', format ='png') 
	else:  
		degree_dist = nx.degree_histogram(G_sub)
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of divergence model - random attack")
		plt.savefig('divergence_duplication_network_%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# For each of the four networks, try to remove vertices in decreasing order of degree. That is, start by removing the vertex with the highest degree, then the vertex with the second highest degree, etc. Plot the degree distribution at regular intervals. What do you observe? How does this type of attach compare to the random attack?
#########################
# Making attacks after sorting the data
#########################
def attacking_sorted(G, n):
	lists = sorted(nx.degree(G).items(), key=lambda (k,v):v, reverse=True)
	for nets in lists[:n]:
		G.remove_nodes_from(nets)
	return G

# random network
print '\nSorted attack in randomized network \n'
for i in range(1, 4):
	G_sub = attacking_sorted(G_random, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_network_sorted' + str(i*100) +'.png', format ='png')  
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of random model - sorted attack")
		plt.savefig('random_network_sorted%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# yeast network
print '\nSorted attack in yeast network \n'
for i in range(1, 4):
	G_sub = attacking_sorted(G_yeast, 200*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig('Graph_yeast_sorted' + str(i*100) +'.png', format ='png')
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of yeast model - sorted attack")
		plt.savefig('yeast_network_sorted%s.png' % (i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# partial duplication
print '\nSorted attack in partial duplication network \n'
for i in range(1, 4):
	G_sub = attacking_sorted(G_partial_dupli, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig("Graph_partial_sorted" + str(i*100) +'.png',  format ='png')
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of partial model - sorted attack")
		plt.savefig('partial_duplication_sorted%s.png' %(i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

# divergence duplication
print '\nSorted attack in divergence duplication network \n'
for i in range(1, 4):
	G_sub = attacking_sorted(G_dupli_div, 100*i)
	if argv[2] == 'images':
		plt.figure()
		nx.draw(G_sub,  node_size = 10)
		plt.savefig("Graph_divergence_sorted" + str(i*100) +'.png', format ='png') 
	else:
		degree_dist = nx.degree_histogram(G_sub) 
		plt.figure()
		plt.loglog(degree_dist,'b-',marker='o')
		plt.title("Degree rank plot of divergence model - sorted attack")
		plt.savefig('div_dupli_network_sorted%s.png' %(i*100))
		print '\nRemoving %d nodes' % (i*100)
		evaluator(G_sub)

#########################
# Making attacks by choosing a specific maximum degree
#########################	
def attacking_degree(G, degree):
	lists = sorted(nx.degree(G_partial_dupli).items(), key=lambda (k,v):v, reverse=True)
	for nets in lists:
		if nets >= degree:
			G.remove_nodes_from(nets)
	return G