import pandas as pd
from shapely import wkt,geometry
import networkx as nx
import numpy as np
import sumolib
import os
import sys 
import multiprocessing 
import random

import graph
sumolib.os.environ["SUMO_HOME"]="/usr/share/sumo"
if 'SUMO_HOME' in sumolib.os.environ:
    tools = sumolib.os.path.join(sumolib.os.environ['SUMO_HOME'], 'tools')
    sumolib.sys.path.append(tools)
else:   
    sumolib.sys.exit("please declare environment variable 'SUMO_HOME'")

EPSILON = 0.0000001

#edgedata = pd.read_csv("tartu_sample/edgedata.csv")
#G = nx.from_pandas_edgelist(edgedata, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)


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
         options.interval_size = str(3600)

    if options.begin_time is None:
         options.begin_time = str(0)

    if options.end_time is None:
         options.begin_time = str(int(options.begin_time) + int(options.number_iteration))

    if options.graphedge_file is None:
        net = sumolib.net.readNet(options.net_file, withInternal=True)
        edgedata = graph.net2graph(net)
        options.graphedge_file = options.input_location + "graphedge.csv"
        edgedata.to_csv(options.graphedge_file, index=False)
    return options

def DAG_all2destination_od(G, destination_id, origins:set):
    def DAG_all2destination(G, destination_id):
        G_reverse = nx.reverse(G)
        ss1 = nx.single_source_dijkstra_path_length(G_reverse,source="d_"+str(destination_id), weight="cost")
        empty =np.inf
        nx.set_node_attributes(G_reverse, empty, "dist_from_destination")
        nx.set_node_attributes(G_reverse, ss1, "dist_from_destination")
    
        def filter_edge(n1, n2):
            if G_reverse.nodes[n1]["dist_from_destination"] ==np.inf:
                return False
            elif G_reverse.nodes[n2]["dist_from_destination"] ==np.inf:
                return False
            else:
                return G_reverse.nodes[n1]["dist_from_destination"] <= G_reverse.nodes[n2]["dist_from_destination"]

        def filter_node(n):
            return  G_reverse.nodes[n]["dist_from_destination"] !=np.inf

        view = nx.subgraph_view(G_reverse, filter_edge=filter_edge, filter_node=filter_node)
        DAG_reverse = view.copy()
        DAG = nx.reverse(DAG_reverse)
        return DAG

    def filter_DAG_origin_based(DAG, origins):
        nodesorted = sorted(DAG.nodes, key=lambda n: DAG.nodes[n]['dist_from_destination'])
        nodesorted.reverse()
        poi = {"o_"+item for item in origins}
        nodedict = dict()
        for node in nodesorted:
            if (node in nodedict.keys()) or (node in poi):
                nodedict[node]=True
                neighbors_temp = list(nx.neighbors(DAG,node))
                for neighbor in neighbors_temp:
                    nodedict[neighbor]=True
            else:
                nodedict[node]=False
                
        def filter_node_DAG(n):
            return  nodedict[n]

        view_DAG = nx.subgraph_view(DAG, filter_node=filter_node_DAG)
        DAG_filtered = view_DAG.copy()
        return DAG_filtered
    ########################################
    if destination_id in origins:
        origins.remove(destination_id)
    DAG = DAG_all2destination(G, destination_id)
    DAG_filtered = filter_DAG_origin_based(DAG, origins)
    sorted_nodes = sorted(DAG_filtered.nodes, key=lambda n: DAG_filtered.nodes[n]['dist_from_destination'])
    return DAG_filtered, sorted_nodes


def calculate_prob_choice(DAG, sorted_nodes, teta, vdict_last):
    def create_vdict(DAG, sorted_nodes):
        neighbors = nx.to_dict_of_lists(DAG)
        neighbors_cost = dict()
        for node in neighbors.keys():
            tempdict = dict()
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
                    print(maxnode)
                    print(vdict[maxnode])
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
            dict[node]["vcost"] = 0
        elif len(nodesset)==1:
            vdict[node]["vcost"] = EPSILON
            mynode = nodesset.pop()
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
                w[vertex] = EPSILON
                if vdict_last!=None:
                    if node in vdict_last.keys():
                        if vertex in vdict_last[node]["neighbors"].keys():
                            w[vertex] = vdict_last[node]["neighbors"][vertex]["prob"]+EPSILON

                vdict[node]["neighbors"][vertex]["prob"] = expdict[vertex]/sumexp
                w[vertex] = w[vertex]*vdict[node]["neighbors"][vertex]["prob"]
                sumw += w[vertex] 

            for vertex in vdict[node]["neighbors"].keys():
                vdict[node]["neighbors"][vertex]["prob"] = w[vertex]/sumw
            vdict[node]["vcost"] = np.abs(-(pi_od*sumexp)/teta)
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

def parallel_route(destination_id, tripdf, G, id2type, option, vdict_last, return_dict):
    origins = set(tripdf[tripdf["to_node"]==destination_id]["from_node"])
    DAG_filtered, sorted_nodes = DAG_all2destination_od(G, destination_id, origins)
    myvdict = None
    if vdict_last!=None:
        myvdict = vdict_last[destination_id]
    dtree = calculate_prob_choice(DAG_filtered, sorted_nodes, teta=float(option.teta), vdict_last=myvdict)
    mydf = decisiontree2path_routedf(dtree, destination_id, tripdf, G, id2type)
    tempdict = dict()
    tempdict["route"] = mydf
    tempdict["vdict"] = dtree
    return_dict[destination_id] = tempdict
def route_df2xml(route, address):
    def convert_row(row):
        return """
        <vehicle id="%s" depart="%s" departLane="%s" departSpeed="%s">
            <route edges="%s"/>
        </vehicle>""" % (row.vehicle_id, row.time, "free", "max", row.path)

    text0 = """<?xml version="1.0" encoding="UTF-8"?>\n\n\n"""
    text1 = """<routes xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/routes_file.xsd">"""
    text2 = '\n'.join(route.apply(convert_row, axis=1))
    text3 = """\n</routes>"""
    with open(address, 'w') as myfile: 
        myfile.write(text0+text1+text2+text3)

def main():
    option = get_options()
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    trip = pd.read_csv(option.trip_file)
    trip = trip[trip["from_node"]!=trip["to_node"]]
    net = sumolib.net.readNet(option.net_file, withInternal=True)
    graphedge = pd.read_csv(option.graphedge_file)
    vdict_last = None
    G = nx.from_pandas_edgelist(graphedge, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)
    id2type = graphedge.set_index("id")["type"].copy()
    processes = []
    destinations = set(trip["to_node"].unique())

    for destination_id in destinations:
        #print(destination_id)
        #parallel_route(destination_id, trip, G, id2type, option, return_dict)
        p = multiprocessing.Process(target = parallel_route, args=(destination_id,trip,G,id2type, option, vdict_last, return_dict))
        p.start()
        processes.append(p)
    for p in processes:
        p.join() 


    route = pd.concat([return_dict[key]["route"] for key in return_dict.keys()])
    vdict = {key:return_dict[key]["vdict"] for key in return_dict.keys()}
    pd.to_pickle(vdict, option.output_location +"vdict.pkl")
    route = route.sort_values(["time"])
    route.to_csv(option.output_location + "route_0.csv", index=False)
    route_df2xml(route, option.output_location + "route_0.xml")
if __name__=="__main__":
    main()