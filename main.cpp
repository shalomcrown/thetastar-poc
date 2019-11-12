#include <iostream>
#include <list>
#include <geos_c.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/IntersectionMatrix.h>
#include <geos/io/WKBReader.h>
#include <geos/io/WKBWriter.h>
#include <geos/io/WKTWriter.h>
#include <geos/util/GeometricShapeFactory.h>
#include <geos/geom/util/SineStarFactory.h>
#include <geos/util/GEOSException.h>
#include <geos/util/IllegalArgumentException.h>
#include <geos/operation/linemerge/LineMerger.h>
#include <geos/operation/polygonize/Polygonizer.h>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include <limits>
#include <vector>
#include <list>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <streambuf>
#include <string>

// proj-data proj-bin libproj-dev

using namespace std;
using namespace geos;
using namespace geos::geom;
using namespace geos::operation::polygonize;
using namespace geos::operation::linemerge;
using geos::util::GEOSException;
using geos::util::IllegalArgumentException;

GeometryFactory::Ptr global_factory;


double zone3points[][2] = {{35.2594967539393,32.1872363481983},{35.2739069684505,32.1869022628542},{35.2741742367258,32.1812450843598},{35.2703879361587,32.1780378650559},{35.2673143509925,32.1782605886187},{35.2666684526605,32.17723606023},{35.2675816192678,32.1759442635659},{35.2706774767903,32.177080153736},{35.271151095,32.175571259},{35.2672475339237,32.1710777537194},{35.2633053268627,32.1700086806181},{35.2584944979068,32.1718127414765},{35.2517014292424,32.1780601374122},{35.2545522908458,32.1867018116477},{35.2594967539393,32.1872363481983}};
double startPointA[] = {35.27083396911621,32.175430856310406};
double endPointA[] = {35.26203632354736, 32.17993497969751};
double snehYaakov[][2] = {{35.2639957699073,32.1808441819468},{35.2656439242718,32.1774142390802},{35.2644412170328,32.1747861010395},{35.2610558188787,32.1749197351772},{35.2583831361255,32.1781269544811},{35.2587394938259,32.1804878242464},{35.2619021684173,32.1790623934447},{35.2619021684173,32.1790623934447},{35.2639957699073,32.1808441819468}};

// --------------------------------------------------------

void wkt_print_geoms(vector<const Geometry*>* geoms)
{
    io::WKTWriter* wkt = new io::WKTWriter();
    for(unsigned int i = 0; i < geoms->size(); i++) {
        const Geometry* g = (*geoms)[i];
        string tmp = wkt->write(g);
        cout << "[" << i << "] (WKT) " << tmp << endl;
    }
    delete wkt;
}

// --------------------------------------------------------

void wkt_print_geoms(list<Geometry *>& geoms)
{
    io::WKTWriter wkt;
    int index = 0;

    for(auto it  = geoms.begin(); it != geoms.end(); it++) {
        string tmp = wkt.write(*it);
        cout << "[" << index++ << "] (WKT) " << tmp << endl;
    }
}
// --------------------------------------------------------

void wkt_print_geom(string name, const Geometry* geom)
{
    io::WKTWriter wkt;

    string tmp = wkt.write(geom);
    cout << "[" << name << "] (WKT) " << tmp << endl;
}


// --------------------------------------------------------

class Vertex;

class Vertex:Point {
    Vertex* parent = nullptr;
    Vertex* localParent = nullptr;
    double gValue = numeric_limits<double>::max();
    double priorityValue = 0;
    double lowerBound = -numeric_limits<double>::max();
    double upperBound = numeric_limits<double>::max();
    Point* wgs = nullptr;
    int id;
};

// --------------------------------------------------------

Polygon* polygonFromArray(double pointArray[][2], int nPoints) {
    CoordinateArraySequence* cl = new CoordinateArraySequence();

    for (int pointIndex = 0; pointIndex < nPoints; pointIndex++) {
        cl->add(Coordinate(pointArray[pointIndex][1], pointArray[pointIndex][0]));
    }

    LinearRing* lr = global_factory->createLinearRing(cl);
    return global_factory->createPolygon(lr, nullptr);
}

// --------------------------------------------------------

string epsgUtmStringFromWgs(Point *point) {
    int epsgNumber = 32630 + (int)ceil(point->getY() / 6.0);

    if (point->getX() < 0) {
        epsgNumber += 100;
    }

    stringstream outout;
    outout << "EPSG:" << epsgNumber;
    return outout.str();
}

int epsgUtmCodeFromWgs(Point *point) {
    int epsgNumber = 32630 + (int)ceil(point->getY() / 6.0);

    if (point->getX() < 0) {
        epsgNumber += 100;
    }

    return epsgNumber;
}

// --------------------------------------------------------

Point *wgsToUtmConv(Point *input, OGRCoordinateTransformation *conversion) {

    double x = input->getX();
    double y = input->getY();

    cout << "Original " << x << "," << y << endl;
    conversion->Transform(1, &y, &x);
    cout << "Transformed " << x << "," << y << endl;

    return global_factory->createPoint(Coordinate(y, x));
}


Polygon *wgsToUtmConv(Polygon *input, OGRCoordinateTransformation *conversion) {

    CoordinateArraySequence* cl = new CoordinateArraySequence();

    for (Coordinate coord : *input->getExteriorRing()->getCoordinates()->toVector()) {

        cout << "Original " << coord.x << "," << coord.y << endl;
        conversion->Transform(1, &coord.y, &coord.x);
        cout << "Transformed " << coord.x << "," << coord.y << endl;

        cl->add(coord);
    }

    LinearRing* lr = global_factory->createLinearRing(cl);
    return global_factory->createPolygon(lr, nullptr);
}

// --------------------------------------------------------

Point *constrainToGrid(Point *point, int gridX, int gridY) {

    int x = point->getX() + gridX / 2;
    int y = point->getY() + gridY / 2;

    x -= (x % gridX);
    y -= (y % gridY);

    return global_factory->createPoint(Coordinate(y, x));
}

// --------------------------------------------------------

list<Vertex> findRoute(list<Polygon *> exclusionZones, Polygon* ffa, Point *startPointWGS, Point *endPointPointWGS, double gridSize) {

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

    list<Polygon *> relevantZones;

    // Filter exclusion zones by intersection
    for (Polygon *zn : exclusionZones) {
        if (zn->intersects(ffa)) {
            relevantZones.push_back(wgsToUtmConv(zn, wgsToUtm));
        }
    }




    CPLFree(wgsToUtm);
    CPLFree(utmToWgs);
    return results;
}

// --------------------------------------------------------

int main(int argc, char **argv) {
    PrecisionModel* pm = new PrecisionModel(PrecisionModel::FLOATING );
    global_factory = GeometryFactory::create(pm, 4326);
    GDALAllRegister();

    list<Geometry *> exclusionZones;

    Polygon* ffa = polygonFromArray(zone3points, sizeof(zone3points) / sizeof(double) / 2);
    Polygon* exclusion1 = polygonFromArray(snehYaakov, sizeof(snehYaakov) / sizeof(double) / 2);

    exclusionZones.push_back(exclusion1);

    wkt_print_geom("FFA", ffa);
    wkt_print_geoms(exclusionZones);

    Point *startPoint = global_factory->createPoint(Coordinate(startPointA[1], startPointA[0]));
    Point *endPointPoint = global_factory->createPoint(Coordinate(endPointA[1], endPointA[0]));

    wkt_print_geom("FFA", ffa);
    wkt_print_geoms(exclusionZones);
    wkt_print_geom("startPoint", startPoint);
    wkt_print_geom("endPointPoint", endPointPoint);

    findRoute(exclusionZones, ffa, startPoint, endPointPoint, 5.0);

    return 0;
}










