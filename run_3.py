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


import graph


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
            if len(neighbors[node])==0 and (sorted_nodes.index(node)!=0):
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
            print("*************************************")
            print(vdict[node])
            print("*************************************")

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

def parallel_route(destination_id, tripdf, G, id2type, option, vdict_last, interval, return_dict):
    origins = set(tripdf[tripdf["to_node"]==destination_id]["from_node"])
    DAG_filtered, sorted_nodes = DAG_all2destination_od(G, destination_id, origins)
    myvdict = None
    if vdict_last!=None:
        myvdict = vdict_last[(interval, destination_id)]
    dtree = calculate_prob_choice(DAG_filtered, sorted_nodes, teta=float(option.teta), vdict_last=myvdict)
    #mydf = decisiontree2path_routedf(dtree, destination_id, tripdf, G, id2type)
    tempdict = dict()
    #tempdict["route"] = mydf
    tempdict["vdict"] = dtree
    return_dict[(interval, destination_id)] = tempdict
    
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

#####  SIMULATION   ####################

def create_edge_data_add(iteration, option):
        
    text0 = """<?xml version="1.0" encoding="UTF-8"?>\n\n\n"""
    text1 = """<additional xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://sumo.dlr.de/xsd/additional_file.xsd">
"""
    text2 = """
        <edgeData id="%s" freq="%s" file="%s" excludeEmpty="%s" withInternal="%s"/>"""%("dump_", option.interval_size, "aggregated_"+str(iteration) +".xml", "true","true")

    text3 = """\n</additional>"""
    with open(option.output_location+"edge_data_add_"+str(iteration) +".xml", 'w') as myfile: 
        myfile.write(text0+text1+text2+text3)

def sumoCommand(iteration, option):
    diroutput = option.output_location
    if (not sumolib.os.path.isdir(diroutput)):
        os.mkdir(diroutput,mode=0o777) 
    sc_dict = dict()
    sc_dict["--net-file"] = option.net_file
    sc_dict["--route-files"] = diroutput + "routes_"+str(iteration) +".xml"
    sc_dict["--statistic-output"] = diroutput + "statistic_output_" + str(iteration) +".xml"
    sc_dict["--additional-files"] = diroutput + "edge_data_add_" + str(iteration) +".xml"
    sc_dict["--tripinfo-output"] = diroutput + "trip_info.xml"
    sc_dict["--tripinfo-output.write-undeparted"] = "true"
    sc_dict["--route-steps"] = "200"
    sc_dict["--no-internal-links"] = "False"
    sc_dict["--time-to-teleport"] = "300"
    sc_dict["--time-to-teleport.highways"] = "0"
    sc_dict["--eager-insert"] = "False"
    sc_dict["--random"] = "False"

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

def cost_update(G_df, net, iteration, option, interval):
    aggregated = pd.read_csv(option.output_location+"aggregated_"+str(iteration) +".csv", sep=";")
    aggregated = aggregated[aggregated["interval_begin"]==interval]
    alledges = set(aggregated["edge_id"])
    aggregated = aggregated[aggregated["edge_traveltime"].notna()]

    _edge2time = aggregated.set_index("edge_id")["edge_traveltime"].to_dict()
    #for key in _edge2time.keys():
    #    edge = net.getEdge(key)
    #    if _edge2time[key] < edge.getLength()/edge.getSpeed():
    #        print(str(key) + "  :  "+str(_edge2time[key]) + "  " + str(edge.getLength()/edge.getSpeed()))
    
    for edge in net.getEdges():
        if not(edge.getID() in _edge2time.keys()):
            _edge2time[edge.getID()] = edge.getLength()/edge.getSpeed()
        else:
            _edge2time[edge.getID()] = max(_edge2time[edge.getID()], edge.getLength()/edge.getSpeed())
            
    _G_df = G_df.copy()
    _G_df["cost"] = G_df.apply(lambda row: _edge2time[row.id] if row.id in _edge2time.keys() else row.cost, axis=1)
    _G_df["speed"] = G_df.apply(lambda row: row.length/row.cost, axis =1)
    return _G_df

def run_simulation(iteration, net, graphedge, option, number_of_interval):
    scommand = sumoCommand(iteration, option)
    create_edge_data_add(iteration, option)
    sumolib.subprocess.call(scommand)
    
    sumolib.subprocess.call(["python3", option.tools+"/xml/xml2csv.py",option.output_location+"aggregated_"+str(iteration) +".xml"])
    G_list = list()

    for interval in range(number_of_interval):
        new_G_df = cost_update(graphedge, net,iteration, option, interval)
        G = nx.from_pandas_edgelist(new_G_df, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)
        G_list.append(G)
    return G_list

def paralle_simulation(iteration, net, trip, graphedge, vdict, number_of_interval, G, id2type, option):
    mydf = list()
    for interval in range(number_of_interval):
        trip_interval = trip[trip["interval"]==interval].reset_index(drop=True)
        destinations = set(trip_interval["to_node"].unique())
            #print(trip_interval)
        for destination_id in destinations:
            mydf.append(decisiontree2path_routedf(vdict[(interval, destination_id)], destination_id, trip_interval, G, id2type))
    route = pd.concat(mydf)
    route = route.sort_values(["time"])
    route = route.reset_index(drop=True)
    route.to_csv(option.output_location + "routes_"+str(iteration) +".xml", index=False)
    route_df2xml(route, option.output_location + "routes_"+str(iteration) +".xml")
    
    scommand = sumoCommand(iteration, option)
    create_edge_data_add(iteration, option)
    sumolib.subprocess.call(scommand)
    
    sumolib.subprocess.call(["python3", option.tools+"/xml/xml2csv.py",option.output_location+"aggregated_"+str(iteration) +".xml"])
    G_list = list()

    for interval in range(number_of_interval):
        new_G_df = cost_update(graphedge, net,iteration, option, interval)
        G = nx.from_pandas_edgelist(new_G_df, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)
        G_list.append(G)
    return G_list



def update_report(iteration, report_fieldnames, option):
    tree = ET.parse(option.output_location+"statistic_output_" + str(iteration)+".xml") 
    root = tree.getroot() 
    report = dict()
    report["iteration"] = iteration
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
    with open(option.output_location+"report.csv", 'a', encoding='UTF8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=report_fieldnames)
        writer.writerow(report)  

def add_interval_for_trips(option, trip):
    interval_size = int(option.interval_size)
    trip["begin_interval"] = trip["time"].apply(lambda x: int(x/interval_size)*interval_size)
    trip["end_interval"] = trip["time"].apply(lambda x: (int(x/interval_size)+1)*interval_size)
    trip["interval"] = trip["time"].apply(lambda x: int(x/interval_size))

    return trip

def main():
    os.system("taskset -p -c 0,1,2,3,4,5,6,7,8,9,10,11,12,13 %d" % os.getpid())
    option = get_options()
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    trip = pd.read_csv(option.trip_file)
    trip = trip[trip["from_node"]!=trip["to_node"]]
    trip = add_interval_for_trips(option, trip)
    #trip.to_csv("tartu-sample/trip_2.csv", index=False)
    net = sumolib.net.readNet(option.net_file, withInternal=True)
    graphedge = pd.read_csv(option.graphedge_file)
    vdict_last = None
    G = nx.from_pandas_edgelist(graphedge, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)
    G_list = list()
    number_of_interval = trip["interval"].max() + 1
    for i in range(number_of_interval):
        G_list.append(G.copy())
    
    id2type = graphedge.set_index("id")["type"].copy()
    report_fieldnames = [ "iteration", "v_loaded", "v_inserted","v_running","v_waiting","t_total",
                     "t_jam","t_yield", "t_wrongLane","s_collisions","s_emergencyStops",
                     "vts_routeLength","vts_speed","vts_duration","vts_waitingTime","vts_timeLoss",
                     "vts_departDelay","vts_departDelayWaiting","vts_totalTravelTime","vts_totalDepartDelay"]
    
    with open(option.output_location +'report.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=report_fieldnames)
        writer.writeheader()

    processes = []

    for iteration in range(int(option.number_iteration)):
        for interval in range(number_of_interval):
            trip_interval = trip[trip["interval"]==interval].reset_index(drop=True)
            destinations = set(trip_interval["to_node"].unique())
            #print(trip_interval)
            for destination_id in destinations:
                #print(destination_id)
                #parallel_route(destination_id, trip, G, id2type, option, return_dict)
                p = multiprocessing.Process(target = parallel_route, args=(destination_id,trip_interval,
                                                    G_list[interval],id2type, option, vdict_last, interval, return_dict))
                p.start()
                processes.append(p)
       
        for p in processes:
            p.join() 

        #route = pd.concat([return_dict[key]["route"] for key in return_dict.keys()])
        vdict_last = {key:return_dict[key]["vdict"] for key in return_dict.keys()}
        #pd.to_pickle(vdict, option.output_location +"vdict.pkl")
        #route = route.sort_values(["time"])
        #route = route.reset_index(drop=True)
        #route.to_csv(option.output_location + "routes_"+str(iteration) +".xml", index=False)
        #route_df2xml(route, option.output_location + "routes_"+str(iteration) +".xml")
        #G_list = run_simulation(iteration, net, graphedge, option, number_of_interval)
        G_list = paralle_simulation(iteration, net, trip, graphedge, vdict_last, number_of_interval, G, id2type, option)
        update_report(iteration, report_fieldnames, option)

if __name__=="__main__":
    main()