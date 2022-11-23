import pandas as pd
from shapely import wkt,geometry
import networkx as nx
import numpy as np
import sumolib
import os
import sys 
import graph
sumolib.os.environ["SUMO_HOME"]="/usr/share/sumo"
if 'SUMO_HOME' in sumolib.os.environ:
    tools = sumolib.os.path.join(sumolib.os.environ['SUMO_HOME'], 'tools')
    sumolib.sys.path.append(tools)
else:   
    sumolib.sys.exit("please declare environment variable 'SUMO_HOME'")


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


def trip2route(trip, id2type):
    temp = list(trip.groupby(["from_node","to_node"]))
    sd = [item[0] for item in temp]
    sd2path = dict()
    for item in sd:
        sd2path[item] = shortest_path(G, source=item[0], target=item[1], id2type=id2type)
    trip["path"] =  trip.apply(lambda row: sd2path[(row.from_node, row.to_node)], axis=1)
    return trip

def main():
    option = get_options()

if __name__=="__main__":
    main()