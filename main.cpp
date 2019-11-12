#include <iostream>
#include <list>
// #include <geos_c.h>
// #include <geos/geom/PrecisionModel.h>
// #include <geos/geom/GeometryFactory.h>
// #include <geos/geom/Geometry.h>
// #include <geos/geom/Point.h>
// #include <geos/geom/LinearRing.h>
// #include <geos/geom/LineString.h>
// #include <geos/geom/Polygon.h>
// #include <geos/geom/GeometryCollection.h>
// #include <geos/geom/Coordinate.h>
// #include <geos/geom/CoordinateSequence.h>
// #include <geos/geom/CoordinateArraySequence.h>
// #include <geos/geom/IntersectionMatrix.h>
// #include <geos/io/WKBReader.h>
// #include <geos/io/WKBWriter.h>
// #include <geos/io/WKTWriter.h>
// #include <geos/util/GeometricShapeFactory.h>
// #include <geos/geom/util/SineStarFactory.h>
// #include <geos/util/GEOSException.h>
// #include <geos/util/IllegalArgumentException.h>
// #include <geos/operation/linemerge/LineMerger.h>
// #include <geos/operation/polygonize/Polygonizer.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Timer.h>

#include "gdal_priv.h"
#include "cpl_conv.h"
#include <limits>
#include <vector>
#include <list>
#include <queue>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <streambuf>
#include <string>
#include <algorithm>

// proj-data proj-bin libproj-dev
// gdal
// libcgal-dev
// libcgal-qt5-dev


using namespace std;
// using namespace geos;
// using namespace geos::geom;
// using namespace geos::operation::polygonize;
// using namespace geos::operation::linemerge;
// using geos::util::GEOSException;
// using geos::util::IllegalArgumentException;
//
// GeometryFactory::Ptr global_factory;


using namespace CGAL;
typedef Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2   Point;
// typedef Polygon_2<Kernel> Polygon;
typedef Kernel::Segment_2 Segment;
typedef Polygon_with_holes_2<Kernel> Polygon;


double zone3points[][2] = {{35.2594967539393,32.1872363481983},{35.2739069684505,32.1869022628542},{35.2741742367258,32.1812450843598},
                            {35.2703879361587,32.1780378650559},{35.2673143509925,32.1782605886187}/*,{35.2666684526605,32.17723606023},
                            {35.2675816192678,32.1759442635659},{35.2706774767903,32.177080153736},{35.271151095,32.175571259},
                            {35.2672475339237,32.1710777537194},{35.2633053268627,32.1700086806181},{35.2584944979068,32.1718127414765},
                            {35.2517014292424,32.1780601374122},{35.2545522908458,32.1867018116477}*//*,{35.2594967539393,32.1872363481983}*/};
double startPointA[] = {35.27083396911621,32.175430856310406};
double endPointA[] = {35.26203632354736, 32.17993497969751};
double snehYaakov[][2] = {{35.2639957699073,32.1808441819468},{35.2656439242718,32.1774142390802},{35.2644412170328,32.1747861010395},
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

Polygon polygonFromArray(double pointArray[][2], int nPoints) {

    Polygon poly;

    for (int pointIndex = 0; pointIndex < nPoints; pointIndex++) {
        poly.outer_boundary().push_back(Point(pointArray[pointIndex][1], pointArray[pointIndex][0]));
    }

    return poly;
}

// --------------------------------------------------------

string epsgUtmStringFromWgs(Point& point) {
    int epsgNumber = 32630 + (int)ceil(to_double(point.x()) / 6.0);

    if (point.y() < 0) {
        epsgNumber += 100;
    }

    stringstream outout;
    outout << "EPSG:" << epsgNumber;
    return outout.str();
}

int epsgUtmCodeFromWgs(Point& point) {
    int epsgNumber = 32630 + (int)ceil(to_double(point.y()) / 6.0);

    if (point.x() < 0) {
        epsgNumber += 100;
    }

    return epsgNumber;
}

// --------------------------------------------------------

Point wgsToUtmConv(Point& input, OGRCoordinateTransformation *conversion) {

    auto x = to_double(input.x());
    auto y = to_double(input.y());

//     cout << "Original " << x << "," << y << endl;
    conversion->Transform(1, &x, &y);
//     cout << "Transformed " << x << "," << y << endl;

    return Point(x, y);
}


Polygon wgsToUtmConv(Polygon& input, OGRCoordinateTransformation *conversion) {

    Polygon output;

    for (auto coord = input.outer_boundary().vertices_begin() ; coord != input.outer_boundary().vertices_end(); coord++) {

        auto x = to_double(coord->x());
        auto y = to_double(coord->y());

//         cout << "Original " << x << "," << y << endl;
        conversion->Transform(1, &x, &y);
//         cout << "Transformed " << x << "," << y << endl;
        output.outer_boundary().push_back(Point(x, y));
    }

    return output;
}

// --------------------------------------------------------

Point constrainToGrid(Point& point, int gridX, int gridY) {

    int x = round(to_double(point.x()) + gridX / 2.0);
    int y = round(to_double(point.y()) + gridY / 2.0);

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

    // Filter exclusion zones by intersection
    for (Polygon zn : exclusionZones) {
        if (do_intersect(zn, ffaWGS)) {
            relevantZones.push_back(wgsToUtmConv(zn, wgsToUtm));
        }
    }


    int vertexCounter = 0;
//     auto comp = [](Vertex& v1, Vertex& v2) {return v1.gValue < v2.gValue;};
//     list<Vertex> openVerticesBacking;
//     priority_queue<Vertex, list<Vertex>, decltype(comp)> openVertices(comp, openVerticesBacking);
//     list<Vertex> closedVertices;

    Vertex start(startPoint, nullptr, vertexCounter++);
    start.parent = &start;
    start.gValue = 0.0;
    start.priorityValue = to_double(squared_distance(startPoint, endPoint));
    cout << "Start: " << start << endl;

    Vertex end(endPoint, nullptr, -1);
    cout << "End: " << end << endl;

//     openVertices.push(start);

//     while (openVertices.top()->gValue < end->gValue /* + h(end), which is zero*/) {
//         Vertex *s = openVertices.top();
//         openVertices.pop();
//         closedVertices.push_back(s);
//
//         Vertex* neighbours[] = {
//             new Vertex(s.x() + gridSize, s.y(), s, vertexCounter++),
//             new Vertex(s.x() - gridSize, s.y(), s, vertexCounter++),
//             new Vertex(s.x(), s.y() + gridSize, s, vertexCounter++),
//             new Vertex(s.x(), s.y() - gridSize, s, vertexCounter++),
//             new Vertex(s.x() + gridSize, s.y() + gridSize, s, vertexCounter++),
//             new Vertex(s.x() + gridSize, s.y() - gridSize, s, vertexCounter++),
//             new Vertex(s.x() - gridSize, s.y() + gridSize, s, vertexCounter++),
//             new Vertex(s.x() - gridSize, s.y() - gridSize, s, vertexCounter++),
//         };
//
//         for (Vertex* neighbour : neighbours) {
//             if (ffa->contains(neighbour->point) == false) {
//                 delete neighbour;
//                 continue;
//             }
//
//             for (Polygon * excl : exclusionZones) {
//                 if (excl->contains(neighbour->point)) {
//                     delete neighbour;
//                     continue;
//                 }
//             }
//
//             cl->clear();
//             global_factory->createLineString(cl);
//
//
//
//             if (find_if(closedVertices.begin(), closedVertices.end(), [neighbour](Vertex *v1) {return *v1 == *neighbour;}) == closedVertices.end()) {
//                 delete neighbour;
//                 continue;
//             }
//
//             if (*neighbour == *end) {
//                 delete neighbour;
//                 neighbour = end;
//             } else {
//                 auto maybeNeighbour = find_if(openVerticesBacking.begin(), openVerticesBacking.end(), [neighbour](Vertex *v1) {return *v1 == *neighbour;});
//
//                 if (maybeNeighbour != openVerticesBacking.end()) {
//                     delete neighbour;
//                     neighbour = *maybeNeighbour;
//                     openVerticesBacking.remove(*maybeNeighbour);
//                 }
//             }
//         }
//     }

    // -----------------------------------
    CPLFree(wgsToUtm);
    CPLFree(utmToWgs);
    return results;
}

// --------------------------------------------------------

int main(int argc, char **argv) {
    GDALAllRegister();

    list<Polygon> exclusionZones;

    Polygon ffa = polygonFromArray(zone3points, sizeof(zone3points) / sizeof(double) / 2);
    Polygon exclusion1 = polygonFromArray(snehYaakov, sizeof(snehYaakov) / sizeof(double) / 2);

    exclusionZones.push_back(exclusion1);

    cout << "FFA: " << ffa << endl;

    Point startPoint(startPointA[1], startPointA[0]);
    Point endPoint(endPointA[1], endPointA[0]);

    cout << "Start:" << startPoint << " End:" << endPoint << endl;

    findRoute(exclusionZones, ffa, startPoint, endPoint, 5);

    return 0;
}










