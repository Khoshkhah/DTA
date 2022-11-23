import pandas as pd
from shapely import wkt,geometry
from shapely.geometry import Point
from shapely.geometry import LineString, Polygon
import networkx as nx
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.environ["SUMO_HOME"], 'tools'))
import sumolib
EPSILON = 0.0000001



def net2edgedf(net):

    def getShape(edge):
        _shape = list()        
        for point in edge.getShape():  
            _shape.append(net.convertXY2LonLat(point[0], point[1]))
        return LineString(_shape)

    def getRawShape(edge):
        _shape = list()        
        for point in edge.getRawShape():
            _shape.append(net.convertXY2LonLat(point[0], point[1]))
        return LineString(_shape)

    edgelist = []
    for edge in net.getEdges():
        box = edge.getBoundingBox()
        xx = (box[0]+box[2])/2
        yy = (box[1]+box[3])/2
        _pt = net.convertXY2LonLat(xx,yy)
        _id = edge.getID()
        if(edge.isSpecial()):
            _from = ""
            _to = ""
        else:
            _from = edge.getFromNode().getID()
            _to = edge.getToNode().getID()

        _type = edge.getType()
        _function = edge.getFunction()
        if(_function!="internal"):
            _function ="real"

        _laneNumber = edge.getLaneNumber()

        _length = edge.getLength()
        _speed = edge.getSpeed()
        _name = edge.getName()
        _shape = getShape(edge)
        _rawshape =getRawShape(edge)
        _bicycle = edge.allows("bicycle")
        _vehicle = edge.allows("passenger")
        _pedestrian = edge.allows("pedestrian")
        _bus = edge.allows("bus")
        edgelist.append({"id":_id, "from":_from, "to":_to, "laneNumber":_laneNumber,
                          "pedestrian_allow":_pedestrian,"vehicle_allow":_vehicle, "bicyle_allow":_bicycle,
                         "bus_allow":_bus,"speed":_speed,"function":_function, "shape":_shape,"rawshape":_rawshape,
                         "length":_length, "name":_name, "type":_type})
        
    edge_df = pd.DataFrame.from_dict(edgelist)
    edge_df["weight"] = edge_df.apply(lambda row: row.length/row.speed, axis=1)
    return edge_df





def net2graph(net):
    def getInteralEdgeid(fromedgeid , toedgeid):
        fromedge = net.getEdge(fromedgeid)
        toedge = net.getEdge(toedgeid)
        conlist = fromedge.getOutgoing()[toedge]
        return net.getLane(conlist[0].getViaLaneID()).getEdge().getID()

    def outgoingInternalEdge(fromedgeid):
        outgoinglist = []
        fromedge = net.getEdge(fromedgeid)
        conlist = fromedge.getOutgoing()
        for edge in conlist:
            toedgeid = edge.getID()
            internaledge = getInteralEdgeid(fromedgeid , toedgeid)
            outgoinglist.append({"fromedge":fromedgeid, "toedge":toedgeid, "internaledge":internaledge})
        return outgoinglist

    edges = net2edgedf(net)
    edges = edges[edges["function"]=="real"]

    #edges = edges[edges["vehicle_allow"]==True]
    nodes = set(edges["from"]).union(set(edges["to"]))
    if "" in nodes:
        nodes.remove("")
        
    realedge=edges[["id"]].copy()
    realedge["from"] = realedge["id"].apply(lambda x: "from_"+x)
    realedge["to"] = realedge["id"].apply(lambda x: "to_"+x)
    realedge["type"] = "real"
    realedge["speed"] = realedge["id"].apply(lambda x: net.getEdge(x).getSpeed())
    realedge["length"] = realedge["id"].apply(lambda x: net.getEdge(x).getLength())
    
    internaledgelist = []
    templist = []
    for fromedgeid in list(edges["id"]):
        templist = outgoingInternalEdge(fromedgeid)
        for item in templist:
            internaledgelist.append(item)   
    internaledge_df = pd.DataFrame(internaledgelist)

    internaledge = internaledge_df.copy()
    internaledge.rename(columns={"fromedge":"from", "toedge":"to", "internaledge":"id"}, inplace=True)
    internaledge["from"] = internaledge["from"].apply(lambda x: "to_"+x)
    internaledge["to"] = internaledge["to"].apply(lambda x: "from_"+x)
    internaledge["type"] = "internal"
    internaledge["speed"] = internaledge["id"].apply(lambda x: net.getEdge(x).getSpeed())
    internaledge["length"] = internaledge["id"].apply(lambda x: net.getEdge(x).getLength())
    
    edge2nodelist = []
    for nodeid in nodes:
        node = net.getNode(nodeid)
        for edge in node.getIncoming():
            if not edge.isSpecial():
                edge2nodelist.append({"fromedge":edge.getID(), "tonode":nodeid})
    edge2node_df = pd.DataFrame(edge2nodelist)

    incomingnode = edge2node_df.copy()
    incomingnode.rename(columns={"fromedge":"from", "tonode":"to"}, inplace=True)
    incomingnode["from"] = incomingnode["from"].apply(lambda x: "to_"+x)
    incomingnode["to"] = incomingnode["to"].apply(lambda x: "d_"+x)
    incomingnode["id"] = incomingnode.apply(lambda row: row["from"] +"_"+row["to"], axis=1)
    incomingnode["type"] = "destination_node"
    incomingnode["speed"] = 1
    incomingnode["length"] = EPSILON
    
    
    node2edgelist = []
    for nodeid in nodes:
        node = net.getNode(nodeid)
        for edge in node.getOutgoing():
            if not edge.isSpecial():
                node2edgelist.append({ "fromnode":nodeid, "toedge":edge.getID()})
    node2edge_df = pd.DataFrame(node2edgelist)

                              
    outgoingnode = node2edge_df.copy()
    outgoingnode.rename(columns={"fromnode":"from", "toedge":"to"}, inplace=True)
    outgoingnode["from"] = outgoingnode["from"].apply(lambda x: "o_"+x)
    outgoingnode["to"]   = outgoingnode["to"].apply(lambda x: "from_"+x)
    outgoingnode["id"]   = outgoingnode.apply(lambda row: row["from"] +"_"+row["to"], axis=1)

    outgoingnode["type"] = "origin_node"
    outgoingnode["speed"] = 1
    outgoingnode["length"] = EPSILON
    
    new_edge_data = pd.concat([realedge, internaledge,incomingnode,outgoingnode])
    new_edge_data["cost"] = new_edge_data["length"]/new_edge_data["speed"]
    
    #G = nx.from_pandas_edgelist(new_edge_data, source="from", target="to", edge_attr=True, create_using=nx.DiGraph)
    return new_edge_data


def get_options(args=None):
    argParser = sumolib.options.ArgumentParser()
    argParser.add_argument("-n", "--net-file", help="the SUMO net filename")
    argParser.add_argument("-o", "--graphoutput-file", help="the graph output filename as edgefile pandas") 

    options = argParser.parse_args()
    #options = argParser.parse_known_args(args=args)[0]

    if options.net_file is None or options.graphoutput_file is None:
        argParser.print_help()
        sys.exit()

    return options


def main():
    option = get_options()
    net = sumolib.net.readNet(option.net_file, withInternal=True)
    edgedata = net2graph(net)
    edgedata.to_csv(option.graphoutput_file, index=False)


if __name__=="__main__":
    main()