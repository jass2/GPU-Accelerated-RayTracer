//
// Created by Jason Seaman on 2019-09-12.
//

#ifndef JASS2_POLYGON_H
#define JASS2_POLYGON_H

#include "Eigen/Dense"
#include <vector>

using namespace std;

typedef Eigen::Matrix <float , 3 , 1 > Vector3d ;

class Polygon {

public:
    int numEdges;
    vector<Vector3d> vertices;
    Vector3d color;
    double kd, ks, shine;
    Vector3d normal;
    double refraction;
    double kt;
    vector<Vector3d> normals;

    Polygon() {
        numEdges = 0;
    }
    Polygon(int edges, vector<Vector3d> vert, Vector3d fill, double kdc, double ksc, double shinec, double refractionc, double t) {
        numEdges = edges;
        vertices = vert;
        color = fill;
        kd = kdc;
        ks = ksc;
        shine = shinec;
        Vector3d edge1 = vertices[0] - vertices[1];
        Vector3d edge2 = vertices[0] - vertices[2];
        normal = edge1.cross(edge2);
        refraction = refractionc;
        kt = t;
        //cout << "c++   " << endl << vertices[0] << ", " << vertices[1] << ", " << vertices[2] << endl;
                //cout << "c++ " << endl << color[0] << ", " << color[1] << ", " << color[2] << endl;

    }
    Polygon(int edges, vector<Vector3d> vert, Vector3d fill, double kdc, double ksc, double shinec, double refractionc, double t, vector<Vector3d> normalsc) {
        numEdges = edges;
        vertices = vert;
        color = fill;
        kd = kdc;
        ks = ksc;
        shine = shinec;
        Vector3d edge1 = vertices[0] - vertices[1];
        Vector3d edge2 = vertices[0] - vertices[2];
        normal = edge1.cross(edge2);
        refraction = refractionc;
        kt = t;
        normals = normalsc;
        //cout << "c++ " << endl << vertices[0] << ", " << vertices[1] << ", " << vertices[2] << endl;
        //cout << "c++ " << endl << color[0] << ", " << color[1] << ", " << color[2] << endl;


    }
    ~Polygon(){}

};




#endif //JASS2_POLYGON_H
