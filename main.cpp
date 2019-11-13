#include <iostream>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Timer.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include "gdal_priv.h"
#include "cpl_conv.h"
#include <limits>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <streambuf>
#include <string>
#include <algorithm>
#include <boost/metaparse/start.hpp>
#include <boost/metaparse/start.hpp>
#include <boost/metaparse/start.hpp>

using namespace std;
using namespace boost::geometry;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

typedef model::d2::point_xy<double> Point;
typedef model::polygon<model::d2::point_xy<double> > Polygon;
typedef model::ring<model::d2::point_xy<double> > Ring;


list<Point> zone3points = {{35.2594967539393,32.1872363481983},{35.2739069684505,32.1869022628542},{35.2741742367258,32.1812450843598},
                            {35.2703879361587,32.1780378650559},{35.2673143509925,32.1782605886187},{35.2666684526605,32.17723606023},
                            {35.2675816192678,32.1759442635659},{35.2706774767903,32.177080153736},{35.271151095,32.175571259},
                            {35.2672475339237,32.1710777537194},{35.2633053268627,32.1700086806181},{35.2584944979068,32.1718127414765},
                            {35.2517014292424,32.1780601374122},{35.2545522908458,32.1867018116477}/*,{35.2594967539393,32.1872363481983}*/};
double startPointA[] = {35.27083396911621,32.175430856310406};
double endPointA[] = {35.26203632354736, 32.17993497969751};
list<Point> snehYaakov = {{35.2639957699073,32.1808441819468},{35.2656439242718,32.1774142390802},{35.2644412170328,32.1747861010395},
                        {35.2610558188787,32.1749197351772},{35.2583831361255,32.1781269544811},{35.2587394938259,32.1804878242464},
                        {35.2619021684173,32.1790623934447},{35.2619021684173,32.1790623934447}/*,{35.2639957699073,32.1808441819468}*/};

// --------------------------------------------------------

// void wkt_print_geoms(vector<const Geometry*>* geoms)
// {
//     io::WKTWriter* wkt = new io::WKTWriter();
//     for(unsigned int i = 0; i < geoms->size(); i++) {
//         const Geometry* g = (*geoms)[i];
//         string tmp = wkt->write(g);
//         cout << "[" << i << "] (WKT) " << tmp << endl;
//     }
//     delete wkt;
// }

// --------------------------------------------------------

// void wkt_print_geoms(list<Polygon *>& geoms)
// {
//     io::WKTWriter wkt;
//     int index = 0;
//
//     for(auto it  = geoms.begin(); it != geoms.end(); it++) {
//         string tmp = wkt.write(*it);
//         cout << "[" << index++ << "] (WKT) " << tmp << endl;
//     }
// }
// --------------------------------------------------------

// void wkt_print_geom(string name, const Geometry* geom)
// {
//     io::WKTWriter wkt;
//
//     string tmp = wkt.write(geom);
//     cout << "[" << name << "] (WKT) " << tmp << endl;
// }


// --------------------------------------------------------

class Vertex;

class Vertex {
public:
    Vertex* parent = nullptr;
    Vertex* localParent = nullptr;
    double gValue = numeric_limits<double>::max();
    double priorityValue = 0;
    double lowerBound = -numeric_limits<double>::max();
    double upperBound = numeric_limits<double>::max();
    bool owned = false;
    Point point;
    int id;

    Vertex(Point pt, Vertex *_parent, int _id) {
        point = pt;
        id = _id;
        parent = _parent;
    }

    Vertex(int x, int y, Vertex *_parent, int _id)  : point(x, y) {
        id = _id;
        parent = _parent;
        owned = true;
    }


    bool operator ==(Vertex &other) {
        return point.x() == other.point.x() && point.y() == other.point.y();
    }

    friend ostream& operator<<(ostream& out, Vertex& d) {
        out << "Vertex:" << d.id << " Parent:" << (d.parent == nullptr ? -1 : d.parent->id) << " gValue:" << d.gValue << " priorityValue:" << d.priorityValue << " (" << d.point.x() << "," << d.point.y()<< ")";
        return out;
    }

    friend ostream& operator<<(ostream& out, Vertex* d) {
        out << "Vertex:" << d->id << " Parent:" << (d->parent == nullptr ? -1 : d->parent->id) << " gValue:" << d->gValue << " priorityValue:" << d->priorityValue << " (" << d->point.x() << "," << d->point.y()<< ")";
        return out;
    }


};

// --------------------------------------------------------

// Polygon polygonFromArray(double pointArray[][2], int nPoints) {
//
//     Polygon poly;
//
//     for (int pointIndex = 0; pointIndex < nPoints; pointIndex++) {
//         poly.boundary().push_back(Point(pointArray[pointIndex][1], pointArray[pointIndex][0]));
//     }
//
//     return poly;
// }

// --------------------------------------------------------

string epsgUtmStringFromWgs(Point& point) {
    int epsgNumber = 32630 + (int)ceil(point.x() / 6.0);

    if (point.y() < 0) {
        epsgNumber += 100;
    }

    stringstream outout;
    outout << "EPSG:" << epsgNumber;
    return outout.str();
}

int epsgUtmCodeFromWgs(Point& point) {
    int epsgNumber = 32630 + (int)ceil(point.y() / 6.0);

    if (point.x() < 0) {
        epsgNumber += 100;
    }

    return epsgNumber;
}

// --------------------------------------------------------

Point wgsToUtmConv(Point& input, OGRCoordinateTransformation *conversion) {

    auto x = input.x();
    auto y = input.y();

//     cout << "Original " << x << "," << y << endl;
    conversion->Transform(1, &x, &y);
//     cout << "Transformed " << x << "," << y << endl;

    return Point(x, y);
}


Polygon wgsToUtmConv(Polygon& input, OGRCoordinateTransformation *conversion) {

    vector<Point> points;

    for (auto coord = input.outer().begin() ; coord != input.outer().end(); coord++) {

        auto x = coord->x();
        auto y = coord->y();

//         cout << "Original " << x << "," << y << endl;
        conversion->Transform(1, &x, &y);
//         cout << "Transformed " << x << "," << y << endl;
        points.push_back(Point(x, y));
    }

    Ring ring(points.begin(), points.end());
    return Polygon({ring});
}

// --------------------------------------------------------

Point constrainToGrid(Point& point, int gridX, int gridY) {

    int x = round(point.x() + gridX / 2.0);
    int y = round(point.y() + gridY / 2.0);

    x -= (x % gridX);
    y -= (y % gridY);

    return Point(x, y);
}

Point constrainToGrid(Point& point, int grid) {
    return constrainToGrid(point, grid, grid);
}

// --------------------------------------------------------

list<Vertex> findRoute(list<Polygon> exclusionZones, Polygon& ffaWGS, Point& startPointWGS, Point& endPointPointWGS, int gridSize) {

    list<Vertex> results;

    OGRSpatialReference wgsSRS; // lat then lon
    wgsSRS.SetWellKnownGeogCS("WGS84");

    OGRSpatialReference utmSRS;
    utmSRS.SetProjCS("UTM");
    utmSRS.importFromEPSG(epsgUtmCodeFromWgs(startPointWGS));

    OGRCoordinateTransformation *wgsToUtm = OGRCreateCoordinateTransformation( &wgsSRS, &utmSRS);
    OGRCoordinateTransformation *utmToWgs = OGRCreateCoordinateTransformation( &utmSRS, &wgsSRS);

    if ( wgsToUtm == nullptr) {
        cout << "Couldn't create transformation " << endl;
    }

    Polygon ffa = wgsToUtmConv(ffaWGS, wgsToUtm);

    Point startPointA =  wgsToUtmConv(startPointWGS, wgsToUtm);
    Point startPoint = constrainToGrid(startPointA, gridSize);

    Point endPointA = wgsToUtmConv(endPointPointWGS, wgsToUtm);
    Point endPoint = constrainToGrid(endPointA, gridSize);

    list<Polygon> relevantZones;

    return results;
    // Filter exclusion zones by intersection
    for (Polygon zn : exclusionZones) {
        if (overlaps(zn, ffaWGS)) {
            relevantZones.push_back(wgsToUtmConv(zn, wgsToUtm));
        }
    }


    int vertexCounter = 0;
    auto comp = [](Vertex& v1, Vertex& v2) {return v1.gValue < v2.gValue;};
    std::set<Vertex, decltype(comp)> openVertices(comp);
    list<Vertex> closedVertices;

    Vertex start(startPoint, nullptr, vertexCounter++);
    start.parent = &start;
    start.gValue = 0.0;
    start.priorityValue = distance(startPoint, endPoint);
    cout << "Start: " << start << endl;

    Vertex end(endPoint, nullptr, -1);
    cout << "End: " << end << endl;

    openVertices.insert(start);

    while (openVertices.begin()->gValue < end.gValue /* + h(end), which is zero*/) {
        Vertex *s = openVertices.top();
        openVertices.pop();
        closedVertices.push_back(s);

        Vertex* neighbours[] = {
            new Vertex(s.x() + gridSize, s.y(), s, vertexCounter++),
            new Vertex(s.x() - gridSize, s.y(), s, vertexCounter++),
            new Vertex(s.x(), s.y() + gridSize, s, vertexCounter++),
            new Vertex(s.x(), s.y() - gridSize, s, vertexCounter++),
            new Vertex(s.x() + gridSize, s.y() + gridSize, s, vertexCounter++),
            new Vertex(s.x() + gridSize, s.y() - gridSize, s, vertexCounter++),
            new Vertex(s.x() - gridSize, s.y() + gridSize, s, vertexCounter++),
            new Vertex(s.x() - gridSize, s.y() - gridSize, s, vertexCounter++),
        };

        for (Vertex* neighbour : neighbours) {
            if (ffa->contains(neighbour->point) == false) {
                delete neighbour;
                continue;
            }

            for (Polygon * excl : exclusionZones) {
                if (excl->contains(neighbour->point)) {
                    delete neighbour;
                    continue;
                }
            }

            cl->clear();
            global_factory->createLineString(cl);



            if (find_if(closedVertices.begin(), closedVertices.end(), [neighbour](Vertex *v1) {return *v1 == *neighbour;}) == closedVertices.end()) {
                delete neighbour;
                continue;
            }

            if (*neighbour == *end) {
                delete neighbour;
                neighbour = end;
            } else {
                auto maybeNeighbour = find_if(openVerticesBacking.begin(), openVerticesBacking.end(), [neighbour](Vertex *v1) {return *v1 == *neighbour;});

                if (maybeNeighbour != openVerticesBacking.end()) {
                    delete neighbour;
                    neighbour = *maybeNeighbour;
                    openVerticesBacking.remove(*maybeNeighbour);
                }
            }
        }
    }

    // -----------------------------------
    CPLFree(wgsToUtm);
    CPLFree(utmToWgs);
    return results;
}

// --------------------------------------------------------

int main(int argc, char **argv) {
    GDALAllRegister();

    list<Polygon> exclusionZones;

    Polygon ffa = Polygon({Ring(zone3points.begin(), zone3points.end())}); //polygonFromArray(zone3points, sizeof(zone3points) / sizeof(double) / 2);
    Polygon exclusion1({Ring(snehYaakov.begin(), snehYaakov.end())});  // polygonFromArray(snehYaakov, sizeof(snehYaakov) / sizeof(double) / 2);

    exclusionZones.push_back(exclusion1);

//     cout << "FFA: " << ffa.outer() << endl;

    Point startPoint(startPointA[1], startPointA[0]);
    Point endPoint(endPointA[1], endPointA[0]);

//     cout << "Start:" << startPoint.outer() << " End:" << endPoint << endl;

    findRoute(exclusionZones, ffa, startPoint, endPoint, 5);

    return 0;
}










