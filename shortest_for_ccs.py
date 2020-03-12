
# Python program to find single source shortest paths
# for Directed Acyclic Graphs Complexity :OV(V+E)
# Reference:
# http://www.utdallas.edu/~sizheng/CS4349.d/l-notes.d/L17.pdf
from step_field_ccs import SteppedFieldCCS
from collections import defaultdict
import numpy as np
from scipy.stats import linregress
import networkx as nx

def get_short_paths(dts_step_fields,
                    pv_step_fields,
                    threshold_n_fields=3,
                    threshold_r2=0.99):
    """
    dts_step_fields: [[drift times in the first field], ..., [drift times in the last field]]
    pv_step_fields: [p/v in the first field, ..., p/v in the last field]
    """

    assert len(dts_step_fields) == len(pv_step_fields), \
        "Please check if len(dts_step_fields) == len(pv_step_fields)"
    
    # if number of fields < threshold_n_fields
    if len(dts_step_fields) < threshold_n_fields:
        print("len(dts_step_fields) < threshold_n_fields:{0}".format(threshold_n_fields))
        return None
    
    num_features = [len(x) for x in dts_step_fields]
    total_num_nodes = sum(num_features)

    # get the node idx for dts_step_fields
    node_idx_step_fields = []
    dts_node_idx = np.zeros(total_num_nodes)
    pv_node_idx = np.zeros(total_num_nodes)
    cum_node_idx = 0
    for i in range(len(dts_step_fields)):
        cum_node_idx_next = cum_node_idx+len(dts_step_fields[i])
        node_idx_step_fields.append([(j + cum_node_idx) for j, x in enumerate(dts_step_fields[i])])
        dts_node_idx[cum_node_idx:cum_node_idx_next] = dts_step_fields[i]
        pv_node_idx[cum_node_idx:cum_node_idx_next] = [pv_step_fields[i] for x in dts_step_fields[i]]
        cum_node_idx = cum_node_idx_next

    # print('node_idx_step_fields',node_idx_step_fields)
    # print('dts_node_idx',dts_node_idx)
    # print('pv_node_idx',pv_node_idx)

    # decide the direction of topological sort
    # to decide the minimum number of start nodes
    direction = 1  ## forward (first --> last)
    iter_list = enumerate(dts_step_fields)
    start_nodes_idx = node_idx_step_fields[0]
    end_nodes_idx = node_idx_step_fields[-1]
    if (len(dts_step_fields[0]) >= len(dts_step_fields[-1])):
        direction = -1  ## backward (last --> first)
        iter_list = reversed(list(enumerate(dts_step_fields)))
        start_nodes_idx = node_idx_step_fields[-1]
        end_nodes_idx = node_idx_step_fields[0]
    
    # define a graph
    g = Graph(total_num_nodes)

    for i in range(1, len(dts_step_fields)):
        add_edges_between_two_fields(g=g,
                                     high_field=dts_step_fields[i-1],
                                     low_field=dts_step_fields[i],
                                     node_idx_for_high_field=node_idx_step_fields[i-1],
                                     node_idx_for_low_field=node_idx_step_fields[i],
                                     high_field_step=pv_step_fields[i-1],
                                     low_field_step=pv_step_fields[i],
                                     direction=direction)
    
    # find shortest paths
    for s in start_nodes_idx:
        print ("Following are shortest distances from source %d " % s)
        path, dist = g.shortestPath(s, end_nodes_idx)
        print(path, dist)
        print(pv_node_idx[path], dts_node_idx[path])
        slope, intercept, r_value, p_value, std_err = linregress(pv_node_idx[path], dts_node_idx[path]/1000)
        print("R2=", r_value**2)

def add_edges_between_two_fields(g,
                                high_field,
                                low_field,
                                node_idx_for_high_field,
                                node_idx_for_low_field,
                                high_field_step,
                                low_field_step,
                                direction=1):
    diff_step = low_field_step - high_field_step
    if direction == 1:
        for i in range(len(high_field)):
            for j in range(len(low_field)):
                if low_field[j] > high_field[i]:
                    # print("g.addEdge:", node_idx_for_high_field[i], node_idx_for_low_field[j],\
                    #     low_field[j], high_field[i], (low_field[j]-high_field[i])/diff_step)
                    g.addEdge(node_idx_for_high_field[i], node_idx_for_low_field[j],\
                        (low_field[j]-high_field[i])/diff_step)
    elif direction == -1:
        for i in range(len(high_field)):
            for j in range(len(low_field)):
                if low_field[j] > high_field[i]:
                    # print("g.addEdge:", node_idx_for_low_field[j], node_idx_for_high_field[i],\
                    #     low_field[j], high_field[i], (low_field[j]-high_field[i])/diff_step)
                    g.addEdge(node_idx_for_low_field[j], node_idx_for_high_field[i],\
                        (low_field[j]-high_field[i])/diff_step)

# Graph is represented using adjacency list. Every
# node of adjacency list contains vertex number of
# the vertex to which edge connects. It contains the drift times.
# Note: a weight of the edge will be decided by start
class Graph:
    def __init__(self,vertices):
 
        self.V = vertices # No. of vertices
 
        # dictionary containing adjacency List
        self.graph = defaultdict(list)
        self.rgraph = defaultdict(list)  # reverse graph for path search
 
    # function to add an edge to graph
    def addEdge(self,u,v,w):
        self.graph[u].append((v,w))
        self.rgraph[v].append(u)
 
 
    # A recursive function used by shortestPath
    # DFS based solution
    def topologicalSortUtil(self,v,visited,stack):
 
        # Mark the current node as visited.
        visited[v] = True
 
        # Recur for all the vertices adjacent to this vertex
        if v in self.graph.keys():
            for node,weight in self.graph[v]:
                if visited[node] == False:
                    self.topologicalSortUtil(node,visited,stack)
 
        # Push current vertex to stack which stores topological sort
        stack.append(v)

    # The function to do Topological Sort. 
    def topologicalSort(self):
         
        # Create a vector to store indegrees of all
        # vertices. Initialize all indegrees as 0.
        in_degree = [0]*(self.V)
         
        # Traverse adjacency lists to fill indegrees of
        # vertices.  This step takes O(V+E) time
        for i in self.graph:
            for j,w in self.graph[i]:
                in_degree[j] += 1
 
        # Create an queue and enqueue all vertices with
        # indegree 0
        queue = []
        for i in range(self.V):
            if in_degree[i] == 0:
                queue.append(i)
 
        #Initialize count of visited vertices
        cnt = 0
 
        # Create a vector to store result (A topological
        # ordering of the vertices)
        top_order = []
 
        # One by one dequeue vertices from queue and enqueue
        # adjacents if indegree of adjacent becomes 0
        while queue:
 
            # Extract front of queue (or perform dequeue)
            # and add it to topological order
            u = queue.pop(0)
            top_order.append(u)
 
            # Iterate through all neighbouring nodes
            # of dequeued node u and decrease their in-degree
            # by 1
            for i,w in self.graph[u]:
                in_degree[i] -= 1
                # If in-degree becomes zero, add it to queue
                if in_degree[i] == 0:
                    queue.append(i)
 
            cnt += 1
 
        # Check if there was a cycle
        if cnt != self.V:
            print ("There exists a cycle in the graph")
        else :
            #Print topological order
            print (top_order)

    def get_shortest_path(self, start, dist, end_nodes):
        min_idx = np.argmin(dist[end_nodes])
        parent = end_nodes[min_idx]
        min_dist = dist[parent]
        path = [parent]
        while 1:
            parents = []
            for p in self.rgraph[parent]:
                parents.append(p)
            if start in parents:
                path.insert(0, start)
                break
            min_idx = np.argmin(dist[parents])
            parent = parents[min_idx]
            path.insert(0, parent)
        return path, min_dist
 
    ''' The function to find shortest paths from given vertex.
        It uses recursive topologicalSortUtil() to get topological
        sorting of given graph.'''
    def shortestPath(self, start, end_nodes):
 
        # Mark all the vertices as not visited
        visited = [False]*self.V
        topological_sort =[]
 
        # Call the recursive helper function to store Topological
        # Sort starting from source vertice
        for i in range(self.V):
            if visited[i] == False:
                self.topologicalSortUtil(start,visited,topological_sort)
        
        # print('stack:', stack)
        print('topological_sort:', topological_sort)
        # print(self.topologicalSort())

        # if no path from a start node to end nodes
        if topological_sort[0] not in end_nodes:
            print('cannot find a valid path.')
            return [],-1

        # Initialize distances to all vertices as infinite and
        # distance to source as 0
        dist = [float("Inf")] * (self.V)
        dist[start] = 0
 
        # Process vertices in topological order
        cur_i= topological_sort.pop()
        pre_i = cur_i
        slope_info = dict()
        # # Update distances from the start
        for node,weight in self.graph[cur_i]:
            # print(i, node,weight)
            slope_info[node] = weight

        while topological_sort:
            # Get the next vertex from topological order
            cur_i = topological_sort.pop()
            # Update distances of all adjacent vertices
            print(cur_i, self.graph[cur_i])
            print('slope_info', slope_info)
            for node,weight in self.graph[cur_i]:
                print(node,weight,dist)
                if pre_i in slope_info:
                    new_dist = abs(slope_info[pre_i] - weight)
                    print(dist[node],dist[pre_i] + new_dist)
                # if cur_i in slope_info:
                #     new_dist = abs(slope_info[cur_i] - weight)
                #     print(dist[node],dist[cur_i] + new_dist)

                    if dist[node] > dist[cur_i] + new_dist:
                        dist[node] = dist[cur_i] + new_dist
                else:
                    # dist[node] = 0
                    dist[cur_i] = 0
                slope_info[node] = weight
            pre_i = cur_i

        # Print the calculated shortest distances
        # for i in range(self.V):
        #     print ("%d: %d" %(i, dist[i])) if dist[i] != float("Inf") else  print ("%d: Inf" %(i))
        dist_arr = np.array(dist)
        print(dist_arr)
        # find the shortest path
        shortest_path, dist = self.get_shortest_path(start, dist_arr, end_nodes)
        return shortest_path, dist


def get_possible_ccs_values(ccs_df,
                            adduct_mass,
                            old_drift_tube_length=90.33,
                            drift_tube_length=90.33,
                            neutral_mass=28.013,
                            threshold_n_fields=3,
                            threshold_r2=0.99):
    
    ccs_df['pv'] = ccs_df.ImsPressure/ccs_df.ImsField

    frames = list(ccs_df.frame.drop_duplicates())
    # print('frames:', frames)
    dts_step_fields = []
    pv_step_fields = []
    ccs_df_idx = []
    for frame in frames[::-1]:
        # print(frame)
        pv = list(ccs_df[ccs_df.frame==frame]['pv'])[0]
        # print('pv:', pv)
        pv_step_fields.append(pv)
        dts_step_fields.append(list(ccs_df[ccs_df.frame==frame]['dt']))
        ccs_df_idx.append(list(ccs_df[ccs_df.frame==frame].index))
    
    # print("dts_step_fields", dts_step_fields)
    # print("pv_step_fields", pv_step_fields)
    # print("ccs_df_idx", ccs_df_idx)

    assert len(dts_step_fields) == len(pv_step_fields), \
        "Please check if len(dts_step_fields) == len(pv_step_fields)"
    
    # if number of fields < threshold_n_fields
    if len(dts_step_fields) < threshold_n_fields:
        print("len(dts_step_fields) < threshold_n_fields:{0}".format(threshold_n_fields))
        return []
    
    num_features = [len(x) for x in dts_step_fields]
    total_num_nodes = sum(num_features)

    # get the node idx for dts_step_fields
    node_idx_step_fields = []
    dts_node_idx = np.zeros(total_num_nodes)
    pv_node_idx = np.zeros(total_num_nodes)
    df_node_idx = np.zeros(total_num_nodes).astype(int)
    cum_node_idx = 0
    for i in range(len(dts_step_fields)):
        cum_node_idx_next = cum_node_idx+len(dts_step_fields[i])
        node_idx_step_fields.append([(j + cum_node_idx) for j, x in enumerate(dts_step_fields[i])])
        dts_node_idx[cum_node_idx:cum_node_idx_next] = dts_step_fields[i]
        pv_node_idx[cum_node_idx:cum_node_idx_next] = [pv_step_fields[i] for x in dts_step_fields[i]]
        df_node_idx[cum_node_idx:cum_node_idx_next] = ccs_df_idx[i]
        cum_node_idx = cum_node_idx_next

    # print(node_idx_step_fields)
    # print(dts_node_idx)
    # print(pv_node_idx)
    # print(df_node_idx)

    # decide the direction of topological sort
    # to decide the minimum number of start nodes
    direction = 1  ## forward (first --> last)
    start_nodes_idx = node_idx_step_fields[0]
    end_nodes_idx = node_idx_step_fields[-1]
    if (len(dts_step_fields[0]) >= len(dts_step_fields[-1])):
        direction = -1  ## backward (last --> first)
        start_nodes_idx = node_idx_step_fields[-1]
        end_nodes_idx = node_idx_step_fields[0]
    # print ("direction",direction, "start_nodes_idx", start_nodes_idx)
    # print ("node_idx_step_fields",node_idx_step_fields)
    # define a graph
    DG = nx.DiGraph()

    weighted_edges = []
    for k in range(1, len(dts_step_fields)):
        diff_step = pv_step_fields[k] - pv_step_fields[k-1]
        diff_step = diff_step**2
        print(k, dts_step_fields[k-1], dts_step_fields[k])
        
        if direction == 1:
            for i in range(len(dts_step_fields[k-1])):
                for j in range(len(dts_step_fields[k])):
                    if dts_step_fields[k][j] > dts_step_fields[k-1][i]:
                        weighted_edges.append((node_idx_step_fields[k-1][i], node_idx_step_fields[k][j],\
                            (dts_step_fields[k][j]-dts_step_fields[k-1][i])**2+diff_step))
        elif direction == -1:
            for i in range(len(dts_step_fields[k-1])):
                for j in range(len(dts_step_fields[k])):
                    if dts_step_fields[k][j] > dts_step_fields[k-1][i]:
                        weighted_edges.append((node_idx_step_fields[k][j],node_idx_step_fields[k-1][i],\
                            (dts_step_fields[k][j]-dts_step_fields[k-1][i])**2+diff_step))
    
    for start, end, length in weighted_edges:
        # You can attach any attributes you want when adding the edge
        DG.add_edge(start, end, length=length)
    # DG.add_weighted_edges_from(weighted_edges)
    ccs_list = []
    for s in start_nodes_idx:
        best_ccs_of_this_start = None
        best_ccs_r2 = 0
        for e in end_nodes_idx:
            try:
                #NOTE: a shortest path doesn't have the best correlation for linear regression
                # path = nx.dijkstra_path(DG,s,e, weight='length')
                # print('shortest:', path)
                allpaths = list(nx.all_simple_paths(DG,s,e))
                # print('all paths:', allpaths)
                ###################### testing
                best_r = 0
                path = allpaths[0]
                for p in allpaths:
                    
                    # leng = 0
                    # for i in range(1, len(p)):
                    #     _cur = p[i]
                    #     _pre = p[i-1]
                    #     leng += (pv_node_idx[_cur]-pv_node_idx[_pre])**2+(dts_node_idx[_cur]-dts_node_idx[_pre])**2
                    
                    slope, intercept, r_value, p_value, std_err = linregress(pv_node_idx[p], dts_node_idx[p])
                    # print(p, 'r2_value:', r_value**2, 'leng:', leng)
                    # print(p, 'r2_value:', r_value**2)
                    if best_r < r_value:
                        best_r = r_value
                        path = p
                print("BEST:", path, 'r2_value:', best_r**2)
                ######################
            except Exception as e:
                continue

            # print ("Following are shortest distances from source %d to end %d " % (s,e))
            # print('path:', path)
            # print(pv_node_idx[path], dts_node_idx[path])
            # print("index of ccs_df", df_node_idx[path])
            if 'num_isotopes' in ccs_df.columns:
                feature_info = ccs_df.loc[df_node_idx[path]][['intensity_org','mass','dt','num_isotopes','mppid','intensity_z','intensity','ImsPressure','ImsTemperature','ImsField','frame']]
            else:
                feature_info = ccs_df.loc[df_node_idx[path]][['intensity_org','mass','dt','mppid','intensity_z','intensity','ImsPressure','ImsTemperature','ImsField','frame']]
            # print(feature_info)
            ccs = SteppedFieldCCS(feature_info, adduct_mass, old_drift_tube_length)
            ccs.compute(drift_tube_length=drift_tube_length, neutral_mass=neutral_mass)
            # print(ccs.to_dict())
            # remove the redundant regression lines which share the same start nodes(features)
            if (ccs.r2 >= threshold_r2) & (best_ccs_r2 < ccs.r2):
                best_ccs_of_this_start = ccs
                best_ccs_r2 = ccs.r2
        if best_ccs_of_this_start != None:
            ccs_list.append(best_ccs_of_this_start)
    return ccs_list

if __name__ == '__main__':
    import pandas as pd
    df = pd.read_csv('test_185.02.txt', sep='\t', index_col=0)
    # df = pd.read_csv('test_163.04.txt', sep='\t', index_col=0)
    get_possible_ccs_values(df, 185.02, old_drift_tube_length=78.12, drift_tube_length=78.236)

    # df = pd.read_csv('test_163.04.txt', sep='\t')
    # df['pv'] = df.ImsPressure/df.ImsField

    # frames = list(df.frame.drop_duplicates())
    # print('frames:', frames)
    # dts_step_fields = []
    # pv_step_fields = []
    # for frame in frames[::-1]:
    #     print(frame)
    #     pv = list(df[df.frame==frame]['pv'])[0]
    #     print('pv:', pv)
    #     pv_step_fields.append(pv)
    #     dts_step_fields.append(list(df[df.frame==frame]['dt']))
    # print(dts_step_fields)
    # print(pv_step_fields)

    # dts_step_fields=[[6.3,10.89172,12.842286,12.074029],
    #                  [14.46759,9.432969],
    #                  [12.668278,13.259455,15.930501,10.739446,19.34601,15.482841],
    #                  [14.332945],
    #                  [17.718655,9.759595,17.06831],
    #                  [13.200224,19.0303]]
    # pv_step_fields=[0.262170969,0.306507365,0.334344946,0.368910295,0.410865282,0.462187688]
    
    # get_short_paths(dts_step_fields, pv_step_fields)