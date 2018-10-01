import pandas as pd
import numpy as np
from shapely.geometry import Point,LineString,Polygon
from shapely.ops import linemerge, unary_union, polygonize, transform
import shapely
import geopandas as gpd
import geopandas.tools
import geopandas
from functools import partial
import pyproj

#READING CELLS FROM CSV FILE
df_cells = pd.read_csv('stations.csv',header=0,sep=";",encoding='utf-8')
df_cells.index = df_cells["CELL"]

df_cells["LATITUDE"] = pd.Series(df_cells["LATITUDE"]).str.replace(",",".")
df_cells["LONGITUDE"]= pd.Series(df_cells["LONGITUDE"]).str.replace(",",".")

#CONVERTING XSL FLOAT FORMAT TO PYTHON FLOATS
cols_to_convert = ['LATITUDE','LONGITUDE']
for col in cols_to_convert:
    df_cells[col] = pd.to_numeric(df_cells[col],errors='coerce')

    
stations = df_cells.copy() 
cells_2 = df_cells.copy() #for first points
cells_3 = df_cells.copy() #for second points
stations_impact_area = df_cells.copy() #for impact area for each station

#CREATING GEOMETRY COLUMN AND CONVERTING LAT LNG TO POINTS
geometry = [Point(xy) for xy in zip(df_cells['LONGITUDE'],df_cells['LATITUDE'])] 
crs = {'init': 'epsg:4326'} #Coordinate system = WGS84
stations = gpd.GeoDataFrame(stations, crs=crs, geometry=geometry) #Creating geodataframe

#CREATING BUFFER FOR EACH CELL POINT
stations_buffer = stations.copy()
stations_buffer["geometry"] = stations_buffer.geometry.buffer(0.1)

#FINDING LAT LON VALUES OF NEW POINTS FOR CALCULATE IMPACT AREA OF EACH CELL
R = 6378.1 #Radius of the earth
d = 500 #Distance in km

    
cells_2["lat1"]= np.radians(cells_2["LATITUDE"])#lat converted to radians
cells_2["lon1"] = np.radians(cells_2["LONGITUDE"]) #lon converted to radians
cells_2["brng"]= cells_2["AZIMUTH"]*np.pi/180 #converted to rad

#NEW LAT/LNG FROM EXISTING CELL POINT
cells_2["lat2"] = np.arcsin(np.sin(cells_2["lat1"])*np.cos(d/R) + np.cos(cells_2["lat1"])*np.sin(d/R)*np.cos(cells_2["brng"]))
cells_2["lon2"] = cells_2["lon1"] + np.arctan2(np.sin(cells_2["brng"])*np.sin(d/R)*np.cos(cells_2["lat1"]),
         np.cos(d/R)-np.sin(cells_2["lat1"])*np.sin(cells_2["lat2"]))

cells_2["lat2"] = np.degrees(cells_2["lat2"])
cells_2["lon2"] = np.degrees(cells_2["lon2"])

brng = 0*np.pi/180
cells_3["lat3"] = np.arcsin(np.sin(cells_2["lat1"])*np.cos(d/R) + np.cos(cells_2["lat1"])*np.sin(d/R)*np.cos(brng))
cells_3["lon3"] = cells_2["lon1"] + np.arctan2(np.sin(brng)*np.sin(d/R)*np.cos(cells_2["lat1"]),
         np.cos(d/R)-np.sin(cells_2["lat1"])*np.sin(cells_3["lat3"]))

cells_3["lat3"] = np.degrees(cells_3["lat3"])
cells_3["lon3"] = np.degrees(cells_3["lon3"])

point2 = [Point(xy) for xy in zip(cells_2["lon2"],cells_2["lat2"])] 
cells_2 = gpd.GeoDataFrame(cells_2, crs=crs, geometry=point2) #Creating geodataframe
point3 = [Point(xy) for xy in zip(cells_3["lon3"],cells_3["lat3"])] 
cells_3 = gpd.GeoDataFrame(cells_3, crs=crs, geometry=point3) #Creating geodataframe

#FUNCTIONS FOR SPATIAL ANALYSES
def create_linestring(point1,point2):
    line = LineString([point1, point2])
    return line

def polygon_area_calculator(polygon):
    geom = polygon
    geom_area = transform(partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                        proj='aea',
                        lat1=geom.bounds[1],
                        lat2=geom.bounds[3])),geom)
    area = geom_area.area        
    return area

def clipped_polygon(azimuth,polygon,line1,line2):
    merged_line = linemerge([line1, line2])
    merged = linemerge([polygon.boundary, merged_line])
    borders = unary_union(merged)
    polygons = polygonize(borders)
    polygon_list = list(polygons)
    if len(polygon_list)>1:
        polygon1 = polygon_list[0]
        polygon2 = polygon_list[1]
        polygon1_area = polygon_area_calculator(polygon=polygon1)
        polygon2_area=polygon_area_calculator(polygon=polygon2)
        polygon_areas = (polygon1_area, polygon2_area)
        polygon_dict = {polygon1_area:polygon1, polygon2_area:polygon2}
        if azimuth >= 180:
            result = polygon_dict[max(polygon_areas)]
        elif azimuth == 0 or azimuth == 360:
            result = polygon_dict[max(polygon_areas)]
        else:
            result = polygon_dict[min(polygon_areas)]
        return result
    else:
        return polygon


#CALCULATION OF IMPACT AREA FOR EACH STATION POINT

polygon_list = []
index_list = []
for index in df_cells.index:
    line1 = create_linestring(point1=stations.loc[index]["geometry"],point2=cells_2.loc[index]["geometry"])
    line2 = create_linestring(point1=stations.loc[index]["geometry"], point2 = cells_3.loc[index]["geometry"])
    polygon = clipped_polygon(azimuth=stations.loc[index]["AZIMUTH"],
                              polygon=stations_buffer.loc[index]["geometry"],
                              line1=line1, line2=line2)
    
    polygon_list.append(polygon)
    index_list.append(index)

impacted_area_data = dict(zip(index_list,polygon_list))
s  = pd.Series(polygon_list,index=index_list)
impacted_area_pd = pd.DataFrame({"CELL": s.index,
                                 "GEOMETRY":s.values})
impacted_area_pd.index = impacted_area_pd["CELL"]
stations_impact_area = gpd.GeoDataFrame(impacted_area_pd, crs=crs, geometry=impacted_area_pd["GEOMETRY"])


#READING COMPLAINTS FROM CSV FILE
df_complaints = pd.read_csv('complaints.csv',header=0,sep=";",encoding='utf-8')
df_complaints.index = df_complaints["Complaint"]

df_complaints["LATITUDE"] = pd.Series(df_complaints["LATITUDE"].str.replace(",","."))
df_complaints["LONGITUDE"] = pd.Series(df_complaints["LONGITUDE"].str.replace(",","."))
cols_to_convert = ['LATITUDE','LONGITUDE']
for col in cols_to_convert:
    df_complaints[col] = pd.to_numeric(df_complaints[col],errors='coerce')
geometry = [Point(xy) for xy in zip(df_complaints['LONGITUDE'],df_complaints['LATITUDE'])] 
crs = {'init': 'epsg:4326'} #Coordinate system = WGS84
gdf_complaints = gpd.GeoDataFrame(df_complaints, crs=crs, geometry=geometry) #Creating geodataframe
gdf_complaints.plot()
complaints = [x for x in gdf_complaints.index]

#FUNCTIONS FOR INTERSECTION BETWEEN POINTS AND IMPACT AREAS OF STATIONS
def intersect_complaints(poly,point):
    return poly.contains(point)

def distance_calculator(poly,point):
    centroid_point = poly.centroid
    distance_pts = centroid_point.distance(point)
    return distance_pts

#RESULT
result = {}
for index in index_list:
    complaint_list = []

    complaint_dict = {}
    impact_area_polygon = stations_impact_area.loc[index]["GEOMETRY"]
    
    distance_list = []
    
    i = 0
    while i < len(complaints):
        for complaint in complaints:
            point = gdf_complaints.loc[complaint]["geometry"]
            intersection = intersect_complaints(poly=impact_area_polygon,point=point)
            if intersection == True:
                distance = distance_calculator(poly=impact_area_polygon,point=point)
                distance_list.append(distance)
                complaint_list.append(complaint)
            i = i+1
    result[index] = complaint_list

print(result) #turns a dict which contains cell IDs and complaints. 
        
    
    