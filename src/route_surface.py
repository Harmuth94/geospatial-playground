import geopandas as gpd
import shapely
import osm2geojson
import requests
from shapely.geometry import LineString, Point


def invert_coords(geometry):
     return geometry.map(lambda polygon: shapely.ops.transform(lambda x, y: (y, x), polygon))

def cut(line, distance, lines):
    # Cuts a line in several segments at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [LineString(coords[: i + 1]), LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            lines.append(LineString(coords[:i] + [(cp.x, cp.y)]))
            line = LineString([(cp.x, cp.y)] + coords[i:])
            if line.length > distance:
                cut(line, distance, lines)
            else:
                lines.append(LineString([(cp.x, cp.y)] + coords[i:]))
            return lines




def get_near_osm_ways(geo, n_splits = 4, buffer_size = 10) -> gpd.GeoDataFrame:
    # splitting to avoid inner rings 

    result = cut(geo, geo.length / n_splits, list())
    
    # buffer 10 meters
    meters_to_degress = 0.0001

    buffer = [
        x.buffer(buffer_size * meters_to_degress).simplify(
            buffer_size * meters_to_degress, True
        )
        for x in result
    ]

    geo = None

    for buf in buffer:

        poly = ""
        for (x, y) in buf.exterior.coords:
            poly += f" {y} {x}"

        query = f"""
        [out:json][timeout:25];
        (
        // query part for: “polyline”
        way(poly:"{poly}");
        );
        out body;
        >;
        out skel qt;"""
        r = requests.post(
            "https://overpass.kumi.systems/api/interpreter",
            data={"data": query},
            timeout=1200,
            headers={
                "Accept-Charset": "utf-8;q=0.7,*;q=0.7",
                "From": "philipharmuth@gmail.com",
                "User-Agent": "outforsk.com",
            },
        )

        if geo is None:
            geo = r.json()
        else:
            geo["elements"] += r.json()["elements"]
    
    
    geojson = osm2geojson.json2geojson(geo)
    osm_gdf = gpd.GeoDataFrame.from_features(geojson["features"])
    osm_gdf.drop_duplicates(subset=["id"], inplace=True)
    return osm_gdf

VALID_SURFACES = [
    "unpaved",
    "compacted",
    "fine_gravel",
    "gravel",
    "pebblestone",
    "dirt",
    "earth",
    "ground",
    "woodchips",
    "clay",
    "sand",
    "grass",
]

HIGHWAYS = ["path", "track", "cycleway","bridleway", "service", "footway"]
HIGHWAY_NON_SURFACES = ["asphalt", "paved", "cobblestone", "concrete", "steel", "paving_stones", "asphalt", "sett"]

def filter_offroad(tags):
    if not isinstance(tags,dict):
        return False
    if tags.get("highway") in HIGHWAYS:
        if tags.get("surface") in HIGHWAY_NON_SURFACES:
            return False
        return True
    if tags.get("surface") in VALID_SURFACES:
        return True
    return False


def filter_road(tags):
    if not isinstance(tags,dict):
        return False
    if tags.get("highway"):
        return not filter_offroad(tags)
    return False



def determine_surface(geometry, paved, unpaved):
    # we do incremental buffer to avoid self-intersection
    # going unpaved -> paved -> unpaved etc.

    paved_part = None

    for buffer_size in [0.000005, 0.00001, 0.00002, 0.00005, 0.0001]:
        if paved_part is None:
            unpaved_part = geometry.intersection(unpaved.buffer(buffer_size))
        else:
            unpaved_part = geometry.difference(paved_part).intersection(unpaved.buffer(buffer_size))
        paved_part = geometry.difference(unpaved_part).intersection(paved.buffer(buffer_size))


    unknown_part = geometry.difference(paved_part).difference(unpaved_part)
    output = gpd.GeoDataFrame(geometry=[unpaved_part, paved_part, unknown_part])
    output["surface"] = ["unpaved", "paved", "unknown"]

    return output

def determine_surface2(geometry, paved, unpaved):
    # we do incremental buffer to avoid self-intersection
    # going unpaved -> paved -> unpaved etc.

    paved_part = []
    unpaved_part = []
    unknown_part = []

    # get line segments
    for idx, (x, y) in enumerate(geometry.coords):
        if idx == len(geometry.coords) - 1:
            break
        next_x, next_y = geometry.coords[idx + 1]
        line = LineString([(x, y), (next_x, next_y)])
        line_length = line.length
        paved_part = None

        # until 90% of the line length is covered by either paved or unpaved
        # we do incremental buffer
        for buffer_size in [0.000005, 0.00001, 0.00002, 0.00005, 0.0001]:
            paved_part = line.intersection(paved.buffer(buffer_size))
            if paved_part.length > line_length * 0.9:
                break
            unpaved_part = line.intersection(unpaved.buffer(buffer_size))
            if unpaved_part.length > line_length * 0.9:
                break
        

    unknown_part = geometry.difference(paved_part).difference(unpaved_part)
    output = gpd.GeoDataFrame(geometry=[unpaved_part, paved_part, unknown_part])
    output["surface"] = ["unpaved", "paved", "unknown"]

    return output


# Read in the data
if __name__ == "__main__":
    # read strava API route stream -> geojson output
    gdf = gpd.read_file("~/git/result.json")

    gdf.geometry = invert_coords(gdf.geometry)
    
    osm_gdf = get_near_osm_ways(gdf.geometry[0], n_splits=4, buffer_size=10)
    
    unpaved = osm_gdf[osm_gdf.tags.apply(filter_offroad) & (osm_gdf.type != "Polygon")].geometry.unary_union.intersection(gdf.geometry.buffer(0.0001))
    
    paved = osm_gdf[(osm_gdf.tags.apply(filter_road)) & (osm_gdf.type != "Polygon")].geometry.unary_union.intersection(gdf.geometry.buffer(0.0001))
    
    marked_roads = determine_surface(gdf.geometry[0], paved[0], unpaved[0])
