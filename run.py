import pandas as pd
from shapely import wkt,geometry
import networkx as nx
import numpy as np
import sumolib
import os
import sys 
import multiprocessing 
import random
import csv
import xml.etree.ElementTree as ET 
import logging
import random

import graph


EPSILON = 0.0000001

#edgedata = pd.read_csv("tartu_sample/edgedata.csv")
#G = nx.from_pandas_edgelist(edgedata, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)

def get_logger(log, name: str, log_dir):
    log.setLevel(logging.INFO)
    if (not os.path.isdir(log_dir)):
        os.mkdir(log_dir,mode=0o777) 
    file_handler = logging.FileHandler(log_dir+ f'/{name}.log')
    formatter = logging.Formatter('%(asctime)s : %(message)s')
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)
    return log


log = logging.getLogger("DTA")

def get_options(args=None):
    argParser = sumolib.options.ArgumentParser()
    argParser.add_argument("-n", "--net-file", help="the SUMO net filename")
    argParser.add_argument("-g", "--graphedge-file", help="graph edge file as csv file") 
    argParser.add_argument("-t", "--trip-file", help="trip file as csv file")
    argParser.add_argument("-ni", "--number-iteration", help="number of iteration for each interval, default is 15")
    argParser.add_argument("-is", "--interval-size", help="the size of each interval")
    argParser.add_argument("-b", "--begin-time", help="the beginning time of the simulation")
    argParser.add_argument("-e", "--end-time", help="the ending time of the simulation")

    argParser.add_argument("-il", "--input-location", help="the location of input files")

    argParser.add_argument("-ol", "--output-location", help="the location of output files")
    argParser.add_argument("-teta", "--teta", help="teta for logit, default is 2.5 ")

    argParser.add_argument("-sample-iteration", "--sample-iteration", help="number of sampling for getting travel time average for each edge, default is 5")
    argParser.add_argument("-sumo-home", "--sumo-home", help="SUMO_HOME")
    argParser.add_argument("-sumo-binary", "--sumo-binary", help="sumo binary")


    options = argParser.parse_args()
    #options = argParser.parse_known_args(args=args)[0]

    if options.net_file is None or options.trip_file is None or options.input_location is None or \
       options.output_location is None:
        argParser.print_help()
        sys.exit()

    if not(options.sumo_home is None):
            sumolib.os.environ["SUMO_HOME"]=options.sumo_home 
    if 'SUMO_HOME' in sumolib.os.environ:
        options.tools = sumolib.os.path.join(sumolib.os.environ['SUMO_HOME'], 'tools')
        sumolib.sys.path.append(options.tools)
    else:   
        sumolib.sys.exit("please declare environment variable 'SUMO_HOME'")

    if options.sumo_binary is None:
        options.sumo_binary = "sumo"

    if (not os.path.isdir(options.output_location)):
        os.mkdir(options.output_location,mode=0o777) 

    if (not os.path.isdir(options.input_location)):
        os.mkdir(options.input_location,mode=0o777) 

    if options.teta is None:
        options.teta = str(2.5)
    
    if options.number_iteration is None:
        options.number_iteration = str(15)

    if options.interval_size is None:
         options.interval_size = str(900)

    if options.begin_time is None:
         options.begin_time = str(0)

    if options.end_time is None:
         options.end_time = str(int(options.begin_time) + int(options.interval_size))

    if options.graphedge_file is None:
        net = sumolib.net.readNet(options.net_file, withInternal=True)
        edgedata = graph.net2graph(net)
        options.graphedge_file = options.input_location + "graphedge.csv"
        edgedata.to_csv(options.graphedge_file, index=False)
    return options

def DAG_all2destination_od(G, destination, origins:set):
    def DAG_all2destination(G, destination):
        G_reverse = nx.reverse(G)
        ss1 = nx.single_source_dijkstra_path_length(G_reverse,source=destination, weight="cost")
        empty =np.inf
        nx.set_node_attributes(G_reverse, empty, "dist_from_destination")
        nx.set_node_attributes(G_reverse, ss1, "dist_from_destination")
    
        def filter_edge(n1, n2):
            if G_reverse.nodes[n1]["dist_from_destination"] ==np.inf:
                return False
            elif G_reverse.nodes[n2]["dist_from_destination"] ==np.inf:
                return False
            elif G_reverse.out_degree(n2)==0:
                return True
            #elif n1==destination:
            #    return True
            else:
                return G_reverse.nodes[n1]["dist_from_destination"] <= G_reverse.nodes[n2]["dist_from_destination"]

        def filter_node(n):
            return  G_reverse.nodes[n]["dist_from_destination"] !=np.inf

        view = nx.subgraph_view(G_reverse, filter_edge=filter_edge, filter_node=filter_node)
        DAG_reverse = view.copy()
        DAG = nx.reverse(DAG_reverse)
        return DAG

    def filter_DAG_origin_based(DAG, origins):
        #nodesorted = sorted(DAG.nodes, key=lambda n: DAG.nodes[n]['dist_from_destination'])
        nodesorted = list(nx.topological_sort(DAG))
        #nodesorted.reverse()
        #poi = {"o_"+item for item in origins}
        nodedict = dict()
        for node in nodesorted:
            if (node in nodedict.keys()) or (node in origins):
                nodedict[node]=True
                neighbors_temp = list(nx.neighbors(DAG,node))
                for neighbor in neighbors_temp:
                    if neighbor not in nodedict.keys():
                        nodedict[neighbor]=True
            else:
                nodedict[node]=False
        nodedict[nodesorted[-1]] = True
                
        def filter_node_DAG(n):
            return  nodedict[n]

        view_DAG = nx.subgraph_view(DAG, filter_node=filter_node_DAG)
        DAG_filtered = view_DAG.copy()
        return DAG_filtered
    ########################################
    if destination in origins:
        origins.remove(destination)
    DAG = DAG_all2destination(G, destination)
 
    # vithout filter nodes
    DAG_filtered = DAG #filter_DAG_origin_based(DAG, origins)
    #sorted_nodes = sorted(DAG_filtered.nodes, key=lambda n: DAG_filtered.nodes[n]['dist_from_destination'])
    sorted_nodes = list(nx.topological_sort(DAG_filtered))
    sorted_nodes.reverse()
    return DAG_filtered, sorted_nodes

#---------------------------------------------------------------------------------------
def vdict2dag(node, meetnode, vdict):

    def filter_DAG_origin_basedx(DAG, origins, destination):
        #nodesorted = sorted(DAG.nodes, key=lambda n: DAG.nodes[n]['dist_from_destination'])
        #nodesorted.reverse()
        nodesorted = list(nx.topological_sort(DAG))
        #nodesorted.reverse()

        origin = origins.pop()
        myindex = nodesorted.index(destination)
        nodedict = dict()
        #print(nodesorted)
        removenode1 = nodesorted[myindex+1:]
        removenode2 = nodesorted[0:nodesorted.index(origin)]

        #print(nodesorted)
        nodesorted = nodesorted[nodesorted.index(origin):myindex+1]
        
        for node in nodesorted:
            if (node in nodedict.keys()) or (node == origin):
                nodedict[node]=True
                neighbors_temp = list(nx.neighbors(DAG,node))
                for neighbor in neighbors_temp:
                    nodedict[neighbor]=True
            else:
                nodedict[node]=False
                
        for node in removenode1:
            nodedict[node]=False
       
        for node in removenode2:
            nodedict[node]=False
                
        def filter_node_DAG(n):
            return  nodedict[n]

        view_DAG = nx.subgraph_view(DAG, filter_node=filter_node_DAG)
        DAG_filtered = view_DAG.copy()
        return DAG_filtered
  
    vedge = []
    #print(vdict)
    for key in vdict.keys():
        if(vdict[key]["order"]<vdict[node]["order"] and vdict[key]["order"]>=vdict[meetnode]["order"]):
            vedge.append((key,vdict[key]['vout'], {'cost':vdict[key]['vcost']})) 
            #print(vdict[key]["order"])
    #print(vedge)
    source = node
    destination = meetnode
    #print(vdict[source]["order"])
    #(vdict[destination]["order"])
    for vertex in vdict[node]["neighbors"].keys():
        vedge.append((source, vertex, {"cost":vdict[node]["neighbors"][vertex]['edgecost']}))
    vg = nx.DiGraph()
    vg.add_edges_from(vedge)
    ss1 = nx.single_source_dijkstra_path_length(vg,source=node, weight="cost")
    empty =np.inf
    nx.set_node_attributes(vg, empty, "dist_from_destination")
    nx.set_node_attributes(vg, ss1, "dist_from_destination")
    vg = filter_DAG_origin_basedx(vg, {node}, meetnode)
    vg = nx.reverse(vg)

    return vg

def vdict2dag_v2(node, meetnode, vdict):
    source = node
    destination = meetnode
    vedge = list()
    for key in vdict.keys():
        #print("---------")
        #print(key)
        for vertex in vdict[key]["neighbors"].keys():
            #print((key, vertex, {"cost":-np.log(vdict[key]["neighbors"][vertex]['prob'])}))
            if vdict[key]["neighbors"][vertex]['prob']==0:
                vdict[key]["neighbors"][vertex]['prob'] = EPSILON
                #print("@@@@@@@@@@@")
                #print(node)
                #print(meetnode)
                #print(key)
                #print(vertex)
                #print("@@@@@@@@@@@")
            vedge.append((key, vertex, {"cost":-np.log(vdict[key]["neighbors"][vertex]['prob'])}))
    #print("=============================")
    vg = nx.DiGraph()
    vg.add_edges_from(vedge)
    return vg


    ##################################
    
def calculate_prob_choice_DAG(DAG, sorted_nodes, teta, vdict_last):
    def create_vdict(DAG, sorted_nodes):
        #sorted_nodes = list(nx.topological_sort(DAG))
        #sorted_nodes.reverse()
        neighbors = nx.to_dict_of_lists(DAG)
        neighbors_cost = dict()
        for node in neighbors.keys():
            tempdict = dict()
            if len(neighbors[node])==0 and (sorted_nodes.index(node)!=0):
                print("if len(neighbors[node])==0 and (sorted_nodes.index(node)!=0):")
                print(node)
            for item in neighbors[node]:
                tempdict[item] = {"edgecost":DAG[node][item]["cost"]}
            neighbors_cost[node] = {"order":sorted_nodes.index(node),
                                    "shortestpathcost":DAG.nodes[node]["dist_from_destination"],
                                    "neighbors":tempdict}
        return neighbors_cost

    def calculate_vcost(vdict, node, teta, vdict_last):
        def find_first_meet(vdict, nodesset):
            if len(nodesset)==0:
                print("Error: there is no node in the nodes set")
                return 0
            while len(nodesset)!=1:
                maxorder = -1
                maxnode = ""
                for node in nodesset:
                    temp_order = vdict[node]["order"]
                    if temp_order > maxorder:
                        maxorder = temp_order
                        maxnode = node
                nodesset.remove(maxnode)
                if "vout" not in vdict[maxnode].keys():
                    print("-----------------------------------------------------------")
                    print(maxnode)
                    print(vdict[maxnode])
                    print("-----------------------------------------------------------")

                else:
                    nodesset.add(vdict[maxnode]["vout"])
            return nodesset.pop()

        def vcost_node2meet(vdict,node,meet):
            cost = 0
            nextnode = node
            while nextnode!=meet:

                cost += vdict[nextnode]["vcost"]
                nextnode = vdict[nextnode]["vout"]

            return cost

        nodesset = set(vdict[node]["neighbors"].keys())
        if len(nodesset)==0:
            vdict[node]["vcost"] = 0
            vdict[node]['vout'] = ""
            #print("*************************************")
            #print(vdict[node])
            #print("*************************************")

        elif len(nodesset)==1:
            mynode = nodesset.pop()
            vdict[node]["vcost"] = vdict[node]["neighbors"][mynode]["edgecost"]
            vdict[node]["vout"] = mynode
            vdict[node]["neighbors"][mynode]["prob"] = 1
        else:  
            meetnode = find_first_meet(vdict, nodesset)
            vg = vdict2dag(node, meetnode, vdict)
            #print("vg = ")
            #display(nx.to_pandas_edgelist(vg))
            #DAG_temp, sorted_nodes_temp = DAG_all2destination_od(vg, node, {meetnode})
            DAG_temp = vg.copy()
            sorted_nodes_temp = list(nx.topological_sort(nx.reverse(vg)))
            #sorted_nodes_temp = sorted_nodes_temp.reverse()
            
            vdict_new = calculate_prob_choice(DAG_temp, sorted_nodes_temp, teta=teta, vdict_last=None)
            vg_new = vdict2dag_v2(meetnode, node, vdict_new)
            sh_path=nx.single_source_dijkstra_path_length(vg_new,source=meetnode, weight="cost")
            vdict[node]["vout"] = meetnode
            vdict[node]['vcost'] = vdict_new[meetnode]['vcost']
            w = dict()
            sumw = 0
            for vertex in vdict[node]["neighbors"].keys():
                w[vertex] = 0.5 #EPSILON 
                if vdict_last!=None:
                    if node in vdict_last.keys():
                        if vertex in vdict_last[node]["neighbors"].keys():
                            w[vertex] = vdict_last[node]["neighbors"][vertex]["prob"]  #+ EPSILON
                if vertex==meetnode:
                        vdict[node]["neighbors"][vertex]["prob"] = np.exp(-vg_new.edges[vertex,node]["cost"])
                else:
                        vdict[node]["neighbors"][vertex]["prob"] = np.exp(-sh_path[vertex])*np.exp(-vg_new.edges[vertex,node]["cost"])
                    
                w[vertex] = w[vertex]*vdict[node]["neighbors"][vertex]["prob"]
                sumw += w[vertex]

            for vertex in vdict[node]["neighbors"].keys():
                vdict[node]["neighbors"][vertex]["prob"] = w[vertex]/sumw
        return {node:vdict[node]}


            
    vdict = create_vdict(DAG, sorted_nodes)
    destination = sorted_nodes[0]
    for node in sorted_nodes:
        #if node!=destination:
        #log.info(f"starting vcost for node {node}")
        vdict.update(calculate_vcost(vdict, node, teta, vdict_last))
        #log.info(f"finished vcost for node {node}")

    return vdict    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#---------------------------------------------------------------------------------------

def calculate_prob_choice(DAG, sorted_nodes, teta, vdict_last):
    def create_vdict(DAG, sorted_nodes):
        neighbors = nx.to_dict_of_lists(DAG)
        neighbors_cost = dict()
        for node in neighbors.keys():
            tempdict = dict()
            if len(neighbors[node])==0 and (sorted_nodes.index(node)!=0):
                print("---------- if len(neighbors[node])==0 and (sorted_nodes.index(node)!=0):")
                print(node)
            for item in neighbors[node]:
                tempdict[item] = {"edgecost":DAG[node][item]["cost"]}
            neighbors_cost[node] = {"order":sorted_nodes.index(node),
                                    "shortestpathcost":DAG.nodes[node]["dist_from_destination"],
                                    "neighbors":tempdict}
        return neighbors_cost

    def calculate_vcost(vdict, node, teta, vdict_last):
        def find_first_meet(vdict, nodesset):
            if len(nodesset)==0:
                print("Error: there is no node in the nodes set")
                return 0
            while len(nodesset)!=1:
                maxorder = -1
                maxnode = ""
                for node in nodesset:
                    temp_order = vdict[node]["order"]
                    if temp_order > maxorder:
                        maxorder = temp_order
                        maxnode = node
                nodesset.remove(maxnode)
                if "vout" not in vdict[maxnode].keys():
                    print("-----------------------------------------------------------")
                    print(maxnode)
                    print(vdict[maxnode])
                    print("-----------------------------------------------------------")

                else:
                    nodesset.add(vdict[maxnode]["vout"])
            return nodesset.pop()

        def vcost_node2meet(vdict,node,meet):
            cost = 0
            nextnode = node
            while nextnode!=meet:

                cost += vdict[nextnode]["vcost"]
                nextnode = vdict[nextnode]["vout"]

            return cost

        nodesset = set(vdict[node]["neighbors"].keys())
        if len(nodesset)==0:
            vdict[node]["vcost"] = 0

        elif len(nodesset)==1:
            mynode = nodesset.pop()
            vdict[node]["vcost"] = vdict[node]["neighbors"][mynode]["edgecost"]
            vdict[node]["vout"] = mynode
            vdict[node]["neighbors"][mynode]["prob"] = 1
        else:  
            meetnode = find_first_meet(vdict, nodesset)
            vdict[node]["vout"] = meetnode
            pi_od = vdict[node]["shortestpathcost"]-vdict[meetnode]["shortestpathcost"]
            expdict = dict()
            sumexp = 0
            for vertex in vdict[node]["neighbors"].keys():
                edgecost = vdict[node]["neighbors"][vertex]["edgecost"]


                cost = edgecost + vcost_node2meet(vdict,vertex,meetnode)# + vdict[meetnode]["shortestpathcost"]

                
                expdict[vertex] = np.exp(-teta*cost/pi_od)
                sumexp += expdict[vertex]
            w = dict()
            sumw = 0
            for vertex in vdict[node]["neighbors"].keys():
                w[vertex] = 0.5 #EPSILON 
                if vdict_last!=None:
                    if node in vdict_last.keys():
                        if vertex in vdict_last[node]["neighbors"].keys():
                            w[vertex] = vdict_last[node]["neighbors"][vertex]["prob"]  #+ EPSILON

                vdict[node]["neighbors"][vertex]["prob"] = expdict[vertex]/sumexp
                w[vertex] = w[vertex]*vdict[node]["neighbors"][vertex]["prob"]
                sumw += w[vertex] 

            for vertex in vdict[node]["neighbors"].keys():
                vdict[node]["neighbors"][vertex]["prob"] = w[vertex]/sumw
            vdict[node]["vcost"] = np.abs(-(pi_od*np.log(sumexp))/teta)
        return vdict


            
    vdict = create_vdict(DAG, sorted_nodes)
    destination = sorted_nodes[0]
    for node in sorted_nodes:
        if node!=destination:
            vdict = calculate_vcost(vdict, node, teta, vdict_last)
    return vdict    


def decisiontree2path_routedf(dtree, destination, tripdf, G, id2type):
    def decisiontree2path_origin(dtree, origin, G, id2type):
        def nodepath2edgepath(path, G):
            path_ = path.copy()
            edgelist = list()
            x = path_[0]
            path_.remove(x)
            for node in path_:
                edgelist.append(G[x][node]["id"])
                x = node
            return " ".join(edgelist)

        def fixpath(path, id2type):
            path = path.split(" ")
            for item in path:
                if id2type[item]!="real":
                    path.remove(item)
            path = " ".join(path)  
            return path

        path = list()
        origin = "o_"+ str(origin)
        path.append(origin)
        nextnodesdict = dtree[origin]["neighbors"]
        while len(nextnodesdict)!=0:
            if len(nextnodesdict)==1:
                nextnode = list(nextnodesdict)[0]
            else:
                nextnode = random.choices(list(nextnodesdict),
                                        [ nextnodesdict[item]["prob"] for item in nextnodesdict.keys()])[0]
            path.append(nextnode)
            nextnodesdict = dtree[nextnode]["neighbors"]
        edgepath = nodepath2edgepath(path, G)
        return fixpath(edgepath, id2type)

    myroutes = tripdf[tripdf["to_node"]==destination].copy().reset_index(drop=True)
    for index,row in myroutes.iterrows():
        newpath = decisiontree2path_origin(dtree, row.from_node, G, id2type)
        myroutes.loc[index,"path"] = newpath
    return myroutes

#**********************************************************888
def decisiontree2path_routedf_v2(vdict_last, origin_id, destination_id, depart_time, number_of_interval , G, id2type, option, net):
    def nodepath2edgepath(path, G):
        path_ = path.copy()
        edgelist = list()
        x = path_[0]
        path_.remove(x)
        for node in path_:
            edgelist.append(G[x][node]["id"])
            x = node
        return " ".join(edgelist)

    def fixpath(path, id2type):
        path = path.split(" ")
        for item in path:
            if id2type[item]!="real":
                path.remove(item)
        path = " ".join(path)  
        return path

    def removeloop(route):
        route = route.split(" ")
        nodes = list()
        mydict = dict()
        end = None
        for edgeid in route:
            edge =  net.getEdge(edgeid)
            begin = edge.getFromNode().getID()
            end = edge.getToNode().getID()
            mydict[(begin, end)] = edgeid
            nodes.append(begin)
        nodes.append(end)

        temp = list()
        for node in nodes:
            if not(node in temp):
                temp.append(node)
            else:
                temp = temp[0:temp.index(node)+1]
        
        x = temp[0]
        temp.remove(x)
        edgelist = list()
        for node in temp:
            edgelist.append(mydict[(x, node)])
            x = node
        return " ".join(edgelist)

    path = list()
    origin = "o_"+ str(origin_id)
    time_path = int(depart_time)
    time_path_list = [time_path]
    interval = min(int(time_path_list[-1]/int(option.interval_size)), number_of_interval-1)

    path.append(origin)
    nextnodesdict = vdict_last[(interval, destination_id)][origin]["neighbors"]
    node = origin
    while len(nextnodesdict)!=0:
        if len(nextnodesdict)==1:
            nextnode = list(nextnodesdict)[0]
        else:
            nextnode = random.choices(list(nextnodesdict),
                                    [ nextnodesdict[item]["prob"] for item in nextnodesdict.keys()])[0]
        time_path += vdict_last[(interval, destination_id)][node]["neighbors"][nextnode]["edgecost"]
        
        #  remove below line if you don't want to update interval for a route
        interval = max(interval, min(int(time_path/int(option.interval_size)), number_of_interval-1))
        node = nextnode
        check = False
        while not (node  in vdict_last[(interval, destination_id)].keys()):
            print("**********")
            print(origin)
            print(depart_time)
            print(interval)
            print(node)
            print(time_path)
            print("**********")
            node = path.pop()
            time_path = time_path_list.pop()

        time_path_list.append(time_path)
        path.append(node)
        nextnodesdict = vdict_last[(interval, destination_id)][node]["neighbors"]


    #path = removeloop(path)
    edgepath = nodepath2edgepath(path, G)
    route = fixpath(edgepath, id2type)
    route = removeloop(route)
  
    return route



#************************************************************
def initialDaG(G,tripdf):
    #log.info(f"number of origins = {len(origins)}")
    destinations = set(tripdf["to_node"])
    DAG_dict = dict()
    for destination_id in destinations:
        origins = set(tripdf[tripdf["to_node"]==destination_id]["from_node"])
        origins = {"o_"+item for item in origins}
        temp = dict()
        temp["DAG"], temp["sorted_nodes"] = DAG_all2destination_od(G, "d_"+destination_id, origins)
        DAG_dict[destination_id] = temp.copy()
    return DAG_dict

def updateDAGcost(DAG_dict, G, destination_id):
    sorted_nodes = DAG_dict[destination_id]["sorted_nodes"].copy()
    DAG_filtered = DAG_dict[destination_id]["DAG"].copy()
    _G = nx.to_pandas_edgelist(G)
    edge2cost = _G.set_index(["source", "target"])["cost"].to_dict()
    #empty = np.inf
    #nx.set_edge_attributes(DAG_filtered, empty, "cost")
    nx.set_edge_attributes(DAG_filtered, edge2cost, "cost")
    ss1 = nx.single_source_dijkstra_path_length( nx.reverse(DAG_filtered),source="d_"+destination_id, weight="cost")
    nx.set_node_attributes(DAG_filtered, ss1, "dist_from_destination")
    #DAG_filtered["speed"] = DAG_filtered["length"]/DAG_filtered["cost"]
    #DAG_filtered = nx.from_pandas_edgelist(DAG_filtered, source="source", target="target", edge_attr=True, create_using=nx.DiGraph)
    return DAG_filtered, sorted_nodes

def parallel_route(destination_id, DAG_dict, G, id2type, option, vdict_last, interval, return_dict):
    #origins = set(tripdf[tripdf["to_node"]==destination_id]["from_node"])
    #log.info(f"number of origins = {len(origins)}")
    #origins = {"o_"+item for item in origins}

    DAG_filtered, sorted_nodes = updateDAGcost(DAG_dict, G, destination_id)
    #pd.to_pickle(DAG_filtered, option.output_location +"DAG_filtered_"+str(destination_id)+"_"+str(interval)+".pkl")
    #pd.to_pickle(sorted_nodes, option.output_location +"sorted_nodes_"+str(destination_id)+"_"+str(interval)+".pkl")

    myvdict = None
    if vdict_last!=None:
        if (interval, destination_id) in vdict_last.keys():
            myvdict = vdict_last[(interval, destination_id)]
    dtree = calculate_prob_choice_DAG(DAG_filtered, sorted_nodes, teta=float(option.teta), vdict_last=myvdict)
    #mydf = decisiontree2path_routedf(dtree, destination_id, tripdf, G, id2type)
    tempdict = dict()
    #tempdict["route"] = mydf
    tempdict["vdict"] = dtree
    #tempdict["DAG"] = DAG_filtered

    return_dict[(interval, destination_id)] = tempdict
    del DAG_filtered, sorted_nodes, tempdict, dtree

def route_df2xml(route, address):
    def convert_row(row):
        return """\n    <vehicle id="%s" depart="%s" departLane="%s" departSpeed="%s">
            <route edges="%s"/>
    </vehicle>""" % (row.vehicle_id, row.time, "free", "max", row.path)

    text0 = """<?xml version="1.0" encoding="UTF-8"?>\n\n\n"""
    text1 = """<routes xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/routes_file.xsd">"""
    text2 = ''.join(route.apply(convert_row, axis=1))
    text3 = """\n</routes>"""
    with open(address, 'w') as myfile: 
        myfile.write(text0+text1+text2+text3)
#####  SIMULATION   ####################

def create_edge_data_add(iteration,sample, option):
        
    text0 = """<?xml version="1.0" encoding="UTF-8"?>\n\n\n"""
    text1 = """<additional xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/additional_file.xsd">
"""
    text2 = """
        <edgeData id="%s" freq="%s" file="%s" excludeEmpty="%s" withInternal="%s"/>"""%("dump_", option.interval_size, "aggregated_"+str(iteration) +"_"+str(sample)+".xml", "true","true")

    text3 = """\n</additional>"""
    with open(option.output_location+"edge_data_add_"+str(iteration) +"_"+str(sample)+".xml", 'w') as myfile: 
        myfile.write(text0+text1+text2+text3)

def sumoCommand(iteration, sample, option):
    diroutput = option.output_location
    if (not sumolib.os.path.isdir(diroutput)):
        os.mkdir(diroutput,mode=0o777) 
    sc_dict = dict()
    sc_dict["--net-file"] = option.net_file
    sc_dict["--route-files"] = diroutput + "routes_"+str(iteration) +"_"+str(sample)+".xml"
    sc_dict["--statistic-output"] = diroutput + "statistic_output_" + str(iteration)+"_"+str(sample) +".xml"
    sc_dict["--additional-files"] = diroutput + "edge_data_add_" + str(iteration) +"_"+str(sample)+".xml"
    sc_dict["--tripinfo-output"] = diroutput + "trip_info_"+"_"+str(sample)+".xml"
    sc_dict["--tripinfo-output.write-undeparted"] = "true"
    sc_dict["--route-steps"] = "200"
    sc_dict["--no-internal-links"] = "False"
    sc_dict["--time-to-teleport"] = "300"
    sc_dict["--time-to-teleport.highways"] = "0"
    sc_dict["--eager-insert"] = "False"
    sc_dict["--random"] = "False"
    #sc_dict["--end"] = "3600"

    sc_dict["--no-warnings"] = "true"
    sc_dict["--no-step-log"] = "true"
    #sc_dict["--end"] = option.end_time
    #sc_dict["--begin"] = option.begin_time
    sumoCmd = list()
    sumoCmd.append(option.sumo_binary)
    for key in sc_dict.keys():
        sumoCmd.append(key)
        sumoCmd.append(sc_dict[key])
    return sumoCmd

def cost_update(G_df, net, iteration, option, interval, best_sample):
    aggregated = pd.read_csv(option.output_location+"aggregated_"+str(iteration)+"_"+str(best_sample) +".csv", sep=";")
    temp_agg = aggregated.groupby("edge_id")["edge_traveltime"].min().reset_index()
    #print(temp_agg)
    _temp_agg = temp_agg.set_index("edge_id")["edge_traveltime"].to_dict()
    aggregated = aggregated[aggregated["interval_begin"].apply(lambda x: int(x)==interval*int(option.interval_size))]
    #log.info("agg len for interval " + str(interval) + " = " + str(len(aggregated)))
    alledges = set(aggregated["edge_id"])
    aggregated = aggregated[aggregated["edge_traveltime"].notna()]
    aggregated["wait_cost"] = aggregated.apply(lambda x: x.edge_waitingTime/(x.edge_left+x.edge_departed+1), axis=1)
    aggregated["edge_traveltime"] = aggregated["edge_traveltime"] #+ aggregated["wait_cost"]
    _edge2time = aggregated.set_index("edge_id")["edge_traveltime"].to_dict()
    #for key in _edge2time.keys():
    #    edge = net.getEdge(key)
    #    if _edge2time[key] < edge.getLength()/edge.getSpeed():
    #        print(str(key) + "  :  "+str(_edge2time[key]) + "  " + str(edge.getLength()/edge.getSpeed()))
    
    for edge in net.getEdges():
        if not(edge.getID() in _edge2time.keys()):
            #if edge.getID() in _temp_agg.keys():
            #    _edge2time[edge.getID()] = max(_temp_agg[edge.getID()], edge.getLength()/edge.getSpeed())
            #else:
            _edge2time[edge.getID()] = edge.getLength()/edge.getSpeed()
        else:
            _edge2time[edge.getID()] = max(_edge2time[edge.getID()], edge.getLength()/edge.getSpeed())
    #print(_edge2time)
    _G_df = G_df.copy()
    _G_df["cost"] = G_df.apply(lambda row: _edge2time[row.id] if row.id in _edge2time.keys() else row.cost, axis=1)
    _G_df["speed"] = G_df.apply(lambda row: row.length/row.cost, axis =1)
    return _G_df

def run_simulation(iteration,sample, net, graphedge, option, number_of_interval):
    scommand = sumoCommand(iteration, option)
    create_edge_data_add(iteration, option)
    sumolib.subprocess.call(scommand)
    
    sumolib.subprocess.call(["python3", option.tools+"/xml/xml2csv.py",option.output_location+"aggregated_"+str(iteration)+"_"+str(sample) +".xml"])
    G_list = list()

    for interval in range(number_of_interval):
        new_G_df = cost_update(graphedge, net,iteration, option, interval)
        G = nx.from_pandas_edgelist(new_G_df, source="source", target="target", edge_attr=True, create_using=nx.DiGraph)
        G_list.append(G)
    return G_list

def parallel_simulation(iteration, sample, net, trip, graphedge, vdict, G, id2type, option, return_G_dict):
    route = trip.copy()  
    temp = list(vdict.keys())
    number_of_interval = len({i for (i,j)in temp})
    #log.info(f"number of interval inside paralle simulation = {number_of_interval}")
    route["path"] = route.apply(lambda row: decisiontree2path_routedf_v2(vdict,
                row.from_node, row.to_node, row.time, number_of_interval , G, id2type, option, net), axis =1)
    #log.info(f"assigning route to trip for sample {sample} is finished.")
    #log.info(f"creating route file for sample {sample}.")

    route.to_csv(option.output_location + "routes_"+str(iteration)+"_"+str(sample) +".csv", index=False)
    route_df2xml(route, option.output_location + "routes_"+str(iteration)+"_"+str(sample) +".xml")
    #log.info(f"route file for sample {sample} is finished")
    scommand = sumoCommand(iteration,sample, option)
    create_edge_data_add(iteration, sample, option)
    sumolib.subprocess.call(scommand)
    sumolib.subprocess.call(["python3", option.tools+"/xml/xml2csv.py",option.output_location+"aggregated_"+str(iteration)+"_"+str(sample) +".xml"])
 

def return_G_best_sample(graphedge, net, iteration, best_sample, option):
    G_list = list()
    aggregated = pd.read_csv(option.output_location+"aggregated_"+str(iteration)+"_"+str(best_sample) +".csv", sep=";")
    number_of_interval = len(aggregated["interval_begin"].unique())
    for interval in range(number_of_interval):
        new_G_df = cost_update(graphedge, net,iteration, option, interval, best_sample)
        G = nx.from_pandas_edgelist(new_G_df, source="source", target="target", edge_attr=True, create_using=nx.DiGraph)
        G_list.append(G)
    return G_list

def update_report(iteration, report_fieldnames, option):
    min_jam = 99999
    min_index = -1
    min_time_cost = np.inf
    min_index_cost = -1

    for sample in range(int(option.sample_iteration)):
        tree = ET.parse(option.output_location+"statistic_output_" + str(iteration)+"_"+str(sample) +".xml") 
        root = tree.getroot() 
        report = dict()
        report["iteration"] = iteration
        report["sample"] = sample

        report["v_loaded"]   = int(root[0].attrib["loaded"])
        report["v_inserted"] = int(root[0].attrib["inserted"])
        report["v_running"]  = int(root[0].attrib["running"])
        report["v_waiting"]  = int(root[0].attrib["waiting"])
        report["t_total"] = int(root[1].attrib["total"])
        report["t_jam"] = int(root[1].attrib["jam"])
        report["t_yield"] = int(root[1].attrib["yield"])
        report["t_wrongLane"] = int(root[1].attrib["wrongLane"])
        report["s_collisions"] = int(root[2].attrib["collisions"])
        report["s_emergencyStops"] = int(root[2].attrib["emergencyStops"])
        report["vts_routeLength"] = float(root[4].attrib["routeLength"])
        report["vts_speed"] = float(root[4].attrib["speed"])
        report["vts_duration"] = float(root[4].attrib["duration"])
        report["vts_waitingTime"] = float(root[4].attrib["waitingTime"])
        report["vts_timeLoss"] = float(root[4].attrib["timeLoss"])
        report["vts_departDelay"] = float(root[4].attrib["departDelay"])
        report["vts_departDelayWaiting"] = float(root[4].attrib["departDelayWaiting"])
        report["vts_totalTravelTime"] = float(root[4].attrib["totalTravelTime"])
        report["vts_totalDepartDelay"] = float(root[4].attrib["totalDepartDelay"])

        if min_time_cost > report["vts_totalTravelTime"] + report["vts_totalDepartDelay"]:
            min_time_cost = report["vts_totalTravelTime"] + report["vts_totalDepartDelay"]
            min_index_cost = sample

        if min_jam > report["t_jam"]:
            min_jam = report["t_jam"]
            min_index = sample
        with open(option.output_location+"report.csv", 'a', encoding='UTF8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=report_fieldnames)
            writer.writerow(report)  
    #return min_index
    return min_index_cost

def add_interval_for_trips(option, trip):
    interval_size = int(option.interval_size)
    trip["begin_interval"] = trip["time"].apply(lambda x: int(x/interval_size)*interval_size)
    trip["end_interval"] = trip["time"].apply(lambda x: (int(x/interval_size)+1)*interval_size)
    trip["interval"] = trip["time"].apply(lambda x: int(x/interval_size))

    return trip

def generateRandoms(a,b,n):
    randoms = []
    for x in range(0,n):
        randoms.append(random.uniform(a,b))
    return randoms

def main():
    os.system("taskset -p -c 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 %d" % os.getpid())
    option = get_options()
    get_logger(log, "run", option.output_location)
    log.info('-------------------------------------------------------------------------')
    #command = "python3 run_7.py -c tartu-sample/tartu.cfg"
    argv ="\"python3 " +" ".join(sys.argv[0:])+"\""
    print(argv)
    log.info("For killing the process use this command:")
    log.info("     pkill -f "+argv)
    #log.info(f"{os.system("taskset -p -c 0,1,2,3,4,5 %d" % os.getpid())}")
    G_list_old = None
    manager = multiprocessing.Manager()
    log.info("reading trips, net and graphedge.")
    trip = pd.read_csv(option.trip_file)
    trip = trip[trip["from_node"]!=trip["to_node"]]
    trip = add_interval_for_trips(option, trip)
    #trip.to_csv("tartu-sample/trip_2.csv", index=False)
    net = sumolib.net.readNet(option.net_file, withInternal=True)
    graphedge = pd.read_csv(option.graphedge_file)
    vdict_last = None
    log.info("creating initial graph list for each interval.")
    G = nx.from_pandas_edgelist(graphedge, source="source", target="target", edge_attr=True, create_using=nx.DiGraph)
    G_list = list()
    number_of_interval = trip["interval"].max() + 1
    for i in range(number_of_interval):
        G_list.append(G.copy())
    
    DAG_dict = initialDaG(G,trip)

    id2type = graphedge.set_index("id")["type"].copy()
    report_fieldnames = [ "iteration", "sample","v_loaded", "v_inserted","v_running","v_waiting","t_total",
                     "t_jam","t_yield", "t_wrongLane","s_collisions","s_emergencyStops",
                     "vts_routeLength","vts_speed","vts_duration","vts_waitingTime","vts_timeLoss",
                     "vts_departDelay","vts_departDelayWaiting","vts_totalTravelTime","vts_totalDepartDelay"]
    
    with open(option.output_location +'report.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=report_fieldnames)
        writer.writeheader()

    destinations = set(trip["to_node"].unique())
    #log.info(f"{list(destinations)}")
    dest_list = []
    packet_size = len(destinations) #10
    for i in range(int(np.ceil(len(destinations)/packet_size))):
        dest_list.append(list(destinations)[i*packet_size:min((i+1)*packet_size, len(destinations))])
        #log.info(f"{dest_list[-1]}")
    for iteration in range(int(option.number_iteration)):
        log.info(f"iteration {iteration} ...")
        log.info("running parallel for calculating vdict for every interval and every destination.")
        
        vdict_new = dict()
        number_of_interval = len(G_list)
        log.info(f"number of interval is {number_of_interval}.")

        for interval in range(len(G_list)):
            #log.info(f"number of destination in interval {interval} is {len(destinations)}")
            
            #log.info(f"creating vdict_new for interval {interval}")
            for new_dests in dest_list:
                processes = []
                return_dict = manager.dict()
                for destination_id in new_dests:
                    #vdic_last = None
                    p = multiprocessing.Process(target = parallel_route, args=(destination_id, DAG_dict,
                                                        G_list[interval],id2type, option, vdict_last, interval, return_dict))
                    p.start()
                    processes.append(p)
        
                for p in processes:
                    p.join() 

            #dag_last = {key:return_dict[key]["DAG"] for key in return_dict.keys()}
            #pd.to_pickle(dag_last, option.output_location +"dag_last.pkl")
            #route = pd.concat([return_dict[key]["route"] for key in return_dict.keys()])
                vdict_new.update({key:return_dict[key]["vdict"] for key in return_dict.keys()})
                del return_dict
        pd.to_pickle(vdict_new, option.output_location +"vdict_"+str(iteration)+".pkl")
        #route = route.sort_values(["time"])
        #route = route.reset_index(drop=True)
        #route.to_csv(option.output_location + "routes_"+str(iteration) +".xml", index=False)
        #route_df2xml(route, option.output_location + "routes_"+str(iteration) +".xml")
        #G_list = run_simulation(iteration, net, graphedge, option, number_of_interval)
        vdict_last = vdict_new.copy()
        del vdict_new
        return_G_dict = manager.dict()
        processes = []
        log.info("starting parallel simulation...")
        for sample in range(int(option.sample_iteration)):
            p = multiprocessing.Process(target = parallel_simulation, args=(iteration, sample, net, trip,
                            graphedge, vdict_last, G, id2type, option, return_G_dict))
            p.start()
            processes.append(p)
       
        for p in processes:
            p.join() 

        
        best_sample = update_report(iteration, report_fieldnames, option)
        log.info(f"best sample is {best_sample}.")
        print(iteration)
        print(best_sample)
        print("---------")
        log.info("creating G_list from best sample.")
        G_list = return_G_best_sample(graphedge, net, iteration, best_sample, option)
        
        for seg in range(len(G_list)):
            g0 = G_list[seg]
            dfg0 = nx.to_pandas_edgelist(g0).sort_values("id").reset_index(drop=True)
            #g1 = G.copy()
            dfg1 = graphedge.sort_values("id").reset_index(drop=True)
            #if seg==0:
            #    dfg0.to_csv(option.output_location+"cost_check_"+str(iteration)+".csv", index=False)
            #dfg0["cost"] = (dfg0["cost"] + 49*dfg1["cost"])/(50)
            #dfg0["cost"] = dfg1["cost"]*generateRandoms(1, 2, len(dfg0))
            #G_list[seg] = nx.from_pandas_edgelist(dfg0, source="source", target="target", edge_attr=True, create_using=nx.DiGraph)

        G_list_old = G_list.copy()

        #pd.to_pickle(G_list, option.output_location +"g_list_"+str(iteration)+".pkl")
        #print("-----------------------------------------------------------------")

if __name__=="__main__":
    main()