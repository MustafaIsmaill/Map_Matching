import osmnx as ox
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import time

place_name = 'leganes tecnologico (fase 1) ,leganes, madrid,spain'

graph = ox.graph_from_place(place_name, network_type='drive')
graph_proj = ox.project_graph(graph)

G = graph_proj

gdf = ox.gdf_from_place(place_name)

geometry = gdf['geometry'].iloc[0]
if isinstance(geometry, Polygon):
    geometry = MultiPolygon([geometry])

geometry_cut = ox.quadrat_cut_geometry(geometry, quadrat_width=0.1)

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

	    circle = gps_point.buffer(distance = 200)

	    # possible_matches_index = list(n_sindex.intersection(circle.bounds))
	    # possible_matches = gdfn.iloc[possible_matches_index]
	    # precise_matches = possible_matches[possible_matches.intersects(circle)]
	    # candidate_nodes = list(precise_matches.index)

	    for poly in geometry_cut:
		    # buffer by the <1 micron dist to account for any space lost in the quadrat cutting
		    # otherwise may miss point(s) that lay directly on quadrat line
		    poly = poly.buffer(1e-14).buffer(0)

		    # find approximate matches with r-tree, then precise matches from those approximate ones
		    possible_matches_index = list(n_sindex.intersection(poly.bounds))
		    possible_matches = gdfn.iloc[possible_matches_index]
		    precise_matches = possible_matches[possible_matches.intersects(poly)]
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

for i in range (0,len(res)):

	plt.plot(res[i][1].x, res[i][1].y, '*')

for i in range(0, len(gps_trajectory['xy'])):
    
	x = gps_trajectory['xy'][i][0]
	y = gps_trajectory['xy'][i][1]

	plt.plot(x, y, '.')
    
plt.show()