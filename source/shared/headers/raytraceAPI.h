#include <iostream>
#include <Eigen>
#include "../polygon.h"
#include "../sphere.h"
#include "../intersection.h"
#include "../../../../Another/Project2/Camera.h"
#include "../ray.h"
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <stdio.h>
#include <CL/cl.hpp>

class raytraceAPI
{
private:
    const double EPSILON = 0.001;
    const int DEPTH_CONST = 5;
    const double inf = numeric_limits<double>::infinity();
    typedef Eigen::Matrix <float, 3, 1 > Vector3d;

public:

    double getLightCoeff(Vector3d normal, Vector3d dottedWith);

    bool raytri(Ray* ray, Vector3d aV, Vector3d bV, Vector3d cV, Vector3d& point1, double t0, double t1);

    bool findHit(Ray* ray, vector<Polygon*> polygons);

    int computeTeapotFile(vector<vector<Vector3d > > pixels, int resX, int resY);

    bool findShadow(Ray* light, vector<Polygon*> polygons);

    Vector3d rayColor(Ray* ray, double t0, double t1, vector<Vector3d> lights, vector<Polygon*> polygons);

    int openNFF(string nffFileName,
        Vector3d& backgroundColor,
        vector<Vector3d>* lights,
        vector<Polygon*>* polygons,
        Camera **camera
        );

    Ray* FireRayAtPixel(Camera* camera, double pixel_x, double pixel_y, Vector3d bgColor, vector<Vector3d> lights, vector<Vector3d> polygons);

    vector<Vector3d> colorPixel(Ray* ray, vector<Vector3d> lights, vector<Polygon*> polygons);

    void deleteAllocated(vector<Polygon*> polygons, vector<Sphere*> spheres);

    int printImageFile(vector<vector<Vector3d > > pixels, Camera* camera);

};