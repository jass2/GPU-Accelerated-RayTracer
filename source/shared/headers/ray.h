//
// Created by Jason Seaman on 2019-09-12.
//

#ifndef JASS2_RAY_H
#define JASS2_RAY_H

#include "Eigen/Dense"
#include <cmath>
#include <vector>
#include "intersection.h"

using namespace std;

typedef Eigen::Matrix <float , 3 , 1 > Vector3d ;

class Ray {

public:
    Vector3d origin;
    Vector3d direction;
    int depth;
    vector<Intersection*> intersections;
    Vector3d background;

    Ray() {}
    Ray(Vector3d sourcec, Vector3d directionc, int depthc, vector<Intersection*> intersectionsc, Vector3d backgroundc) {
        origin = sourcec;
        direction = directionc;
        depth = depthc;
        intersections = intersectionsc;
        background = backgroundc;
    }

    Ray(Vector3d sourcec, Vector3d directionc, int depthc, Vector3d backgroundc) {
        origin = sourcec;
        direction = directionc;
        depth = depthc;
        background = backgroundc;
    }

    ~Ray(){}

};




#endif
