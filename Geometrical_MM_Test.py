import osmnx as ox
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import time

P = Point()
Px = 436596.2936668355
Py = 4467755.533093071

Pxy = (Px, Py)

P = Pxy

place_name = 'leganes tecnologico (fase 1) ,leganes, madrid,spain'

graph = ox.graph_from_place(place_name, network_type='drive')
# graph.add_node(5584375242, osmid = '5584375242', x = Px, y = Py, geometry = P)

graph_proj = ox.project_graph(graph)

G = graph_proj
gdf = ox.gdf_from_place(place_name)
print(type(graph_proj))
nodes, edges = ox.graph_to_gdfs(graph_proj, nodes=True, edges=True)

fig, ax = plt.subplots()

def GeoMM(traj, gdfn, gdfe):

	traj = pd.DataFrame(traj, columns=['timestamp', 'xy'])
	traj['geom'] = traj.apply(lambda row: Point(row.xy), axis=1)
	traj = gpd.GeoDataFrame(traj, geometry=traj['geom'], crs='EPSG3740')
	traj.drop('geom', axis=1, inplace=True)

	#creates an R-tree spatial index
	n_sindex = gdfn.sindex

	res = []

	for gps_frame in traj.itertuples():
	    
	    #integer data type
	    time_stamp = gps_frame[1]

	    #point of shapely.point data type use .x and .y to get x and y float values
	    gps_point = gps_frame[3]

	    buffer_size = 200

	    # print("############")
	    while True:
		    circle = gps_point.buffer(distance = buffer_size)
		    possible_matches_index = list(n_sindex.intersection(circle.bounds))
		    possible_matches = gdfn.iloc[possible_matches_index]
		    precise_matches = possible_matches[possible_matches.intersects(circle)]
		    candidate_nodes = list(precise_matches.index)
		    
		    if len(precise_matches) < 2:
		    	buffer_size = buffer_size + 50
		    else:
		    	break

	    x, y = circle.exterior.coords.xy

	    for i in range(0, len(x)):
	    	plt.plot(x[i], y[i], 'g.')
	    
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
	        dist.append([ls.distance(gps_point), edge, ls])

	        d = ls.distance(gps_point)
	        lenn.append(d.iloc[0])

	    try:
	    	_, idx = min( (lenn[i],i) for i in xrange(len(lenn)) )
	    	true_edge = dist[idx][1]
	    	true_edge_geom = dist[idx][2].item()
	    	pp = true_edge_geom.interpolate(true_edge_geom.project(gps_point)) # projected point
	    	res.append((time_stamp, pp, true_edge, true_edge_geom))
	    except:
	    	pass

	return res

gps_trajectory = {'timestamp': [1,2,3,4,5], 'xy': [(436986, 4.46758e+06), (436931, 4.46755e+06),
 (436863, 4.46754e+06), (436754, 4.46757e+06), (436669, 4.46755e+06)]}

start = time.time()

res = GeoMM(gps_trajectory, nodes, edges)

end = time.time()

print("#Time = %f" % (end - start) )

edges.plot(ax=ax)
# nodes.plot(ax=ax)
# print(nodes[:].x[5584375242], nodes[:].y[5584375242], 'r*')

for i in range (0,len(res)):

	plt.plot(res[i][1].x, res[i][1].y, 'k*')
	# print(i)

for i in range(0, len(gps_trajectory['xy'])):
    
	x = gps_trajectory['xy'][i][0]
	y = gps_trajectory['xy'][i][1]

	plt.plot(x, y, 'r.')
    
plt.show()

# print(graph.number_of_nodes())
# print(nodes.crs)

# x_list = []
# y_list = []

# for i in range(0,len(edges)):

# 	x_points = edges[:].geometry[i].xy[0]
# 	y_points = edges[:].geometry[i].xy[1]

# 	for xy in range(0,len(x_points)):

# 		x_list.append(x_points[xy])
# 		y_list.append(y_points[xy])

# x_tup = tuple(x_list)
# y_tup = tuple(y_list)

# gdf_nodes = gpd.GeoDataFrame(data={'x':x_tup, 'y':y_tup})
# gdf_nodes.crs = gdf.crs
# gdf_nodes.name = 'nodes'
# gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point((row['x'], row['y'])), axis=1)

# nan = float('nan')

# attributes = {k:v for k, v in G[u][v][k].items() if k != 'key'}
# G.add_edge(u, v, key=k+1, **attributes)

# edgelist= [(0,1), (1,2), (2,3)]