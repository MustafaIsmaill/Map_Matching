import osmnx as ox
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import time

place_name = 'leganes tecnologico (fase 1) ,leganes, madrid,spain'

graph = ox.graph_from_place(place_name, network_type='drive')
graph_proj = ox.project_graph(graph)

G = graph_proj

# ox.plot_graph(graph_proj)

nodes, edges = ox.graph_to_gdfs(graph_proj, nodes=True, edges=True)

def GeoMM(traj, gdfn, gdfe):

	traj = pd.DataFrame(traj, columns=['timestamp', 'xy'])
	traj['geom'] = traj.apply(lambda row: Point(row.xy), axis=1)
	traj = gpd.GeoDataFrame(traj, geometry=traj['geom'], crs='EPSG3740')
	traj.drop('geom', axis=1, inplace=True)

	n_sindex = gdfn.sindex

	res = []

	for gps in traj.itertuples():
	    tm = gps[1]
	    p = gps[3]
	    circle = p.buffer(200)
	    possible_matches_index = list(n_sindex.intersection(circle.bounds))
	    possible_matches = gdfn.iloc[possible_matches_index]
	    precise_matches = possible_matches[possible_matches.intersects(circle)]
	    candidate_nodes = list(precise_matches.index)

	    candidate_edges = []
	    for nid in candidate_nodes:
	        candidate_edges.append(G.in_edges(nid))
	        candidate_edges.append(G.out_edges(nid))

	    candidate_edges = [item for sublist in candidate_edges for item in sublist]
	    dist = []
	    lenn = []
	    for edge in candidate_edges:
	        # get the geometry
	        ls = gdfe[(gdfe.u == edge[0]) & (gdfe.v == edge[1])].geometry
	        dist.append([ls.distance(p), edge, ls])

	        d = ls.distance(p)
	        lenn.append(d.iloc[0])

	    _, idx = min( (lenn[i],i) for i in xrange(len(lenn)) )

	    true_edge = dist[idx][1]
	    true_edge_geom = dist[idx][2].item()
	    pp = true_edge_geom.interpolate(true_edge_geom.project(p)) # projected point
	    res.append((tm, pp, true_edge, true_edge_geom))

	return res

gps_trajectory = {'timestamp': [1,2,3,4,5], 
'xy': [(436940, 4.46774e+06), (436882, 4.46774e+06), (436800, 4.46781e+06), (436750, 4.46779e+06), (436704, 4.46774e+06)]}

start = time.time()

res = GeoMM(gps_trajectory, nodes, edges)

end = time.time()

print("#Time = %f" % (end - start) )

fig, ax = plt.subplots()
edges.plot(ax=ax)

for i in range (0,len(res)):

	plt.plot(res[i][1].x, res[i][1].y, '*')

for i in range(0, len(gps_trajectory['xy'])):
    
	x = gps_trajectory['xy'][i][0]
	y = gps_trajectory['xy'][i][1]

	plt.plot(x, y, '.')
    
plt.show()