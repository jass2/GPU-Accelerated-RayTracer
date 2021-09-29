//
// Created by Jason Seaman on 2019-09-03.
//

#include "../../../../Another/Project2/raytraceAPI.h"

#define M_PI = 3.14159262

using namespace std;
using namespace cl;

namespace stdfs = filesystem;

// distance formula to organize intersections on distance.
double distance(const Intersection* j) {
    return sqrt((j->point(0) - j->origin(0)) * (j->point(0) - j->origin(0)) + (j->point(1) - j->origin(1)) * (j->point(1) - j->origin(1)) + (j->point(2) - j->origin(2)) * (j->point(2) - j->origin(2)));
}


double raytraceAPI::getLightCoeff(Vector3d normal, Vector3d dottedWith) {
    double num = normal.dot(dottedWith);
    if (num > 0) {
        return num;
    }
    else {
        return 0;
    }
}


    // taken from textbook section 4.3 /// project 1
    bool raytraceAPI::raytri(Ray* ray, Vector3d aV, Vector3d bV, Vector3d cV, Vector3d& point1, double t0, double t1) {

        double a, b, c, d, e, f, g, h, i, j, k, l;

        a = aV(0) - bV(0);
        b = aV(1) - bV(1);
        c = aV(2) - bV(2);
        d = aV(0) - cV(0);
        e = aV(1) - cV(1);
        f = aV(2) - cV(2);
        g = ray->direction(0);
        h = ray->direction(1);
        i = ray->direction(2);
        j = aV(0) - ray->origin(0);
        k = aV(1) - ray->origin(1);
        l = aV(2) - ray->origin(2);


        // calculations for Cramer's rule

        double m = (a * (e * i - h * f)) + (b * (g * f - d * i)) + (c * (d * h - e * g));

        double t = -1 * (((f * (a * k - j * b)) + (e * (j * c - a * l)) + (d * (b * l - k * c))) / m);
        if (t < t0 || t > t1) {
            return false;
        }


        double r = ((i * (a * k - j * b)) + (h * (j * c - a * l)) + (g * (b * l - k * c))) / m;
        if (r < 0 || r > 1) {
            return false;
        }


        double bee = ((j * (e * i - h * f) + (k * (g * f - d * i)) + (l * (d * h - e * g))) / m);
        if (bee < 0 || (bee > 1 - r)) {
            return false;
        }

        point1 = aV + bee * (bV - aV) + r * (cV - aV);

        return true;
    }

    // new function to test if a ray intersects with any polygon
    bool raytraceAPI::findHit(Ray* ray, vector<Polygon*> polygons) {
        bool intersect = false;
        for (int x = 0; x < polygons.size(); x++) {
            Polygon* p = polygons[x];
            vector <Vector3d> vertices = p->vertices;
            Vector3d point(0, 0, 0);

            for (int y = 0; y < vertices.size() - 2; y++) {
                // calls the function from project one
                intersect = raytri(ray, vertices[0], vertices[y + 1], vertices[y + 2], point, 0, inf);
                if (intersect) {
                    // creates new intersection object (hitrecord)
                    ray->intersections.push_back(new Intersection(point, p, ray->origin));
                    intersect = false;
                }
            }
        }
        // organizes intersections on distance via lambdasont


        if (ray->intersections.size() > 0) {
            //stackoverflow.com/questions/9706517/sort-a-vector-of-objects-by-an-objects-attribute
            sort(ray->intersections.begin(), ray->intersections.end(), [](const Intersection* i, const Intersection* j) {
                return distance(i) < distance(j);
                });
        }
        return ray->intersections.size() > 0;
    }


    int raytraceAPI::computeTeapotFile(vector<vector<Vector3d > > pixels, int resX, int resY) {
        ofstream myfile;
        myfile.open("teapot.ppm");
        myfile << "P3\n" << resX << " " << resY << "\n" << 255 << endl;
        cout << "Building file..." << endl;
        for (int i = 0; i < resY; i++) {
            for (int j = 0; j < resX; j++) {
                for (int x = 0; x < 3; x++) {
                    myfile << (int)pixels[j][i](x) << endl;
                }
            }
        }
        myfile.close();
        return 0;
    }


    // test for first intersection between point of intersection and light source
    bool raytraceAPI::findShadow(Ray* light, vector<Polygon*> polygons) {
        for (int i = 0; i < polygons.size(); i++) {
            bool intersect = false;
            Polygon* p = polygons[i];
            vector <Vector3d> vertices = p->vertices;
            Vector3d point(EPSILON, EPSILON, EPSILON);

            for (int y = 0; y < vertices.size() - 2; y++) {
                intersect = raytri(light, vertices[0], vertices[y + 1], vertices[y + 2], point, 0, inf);
                if (intersect) {
                    return true;
                }
            }
        }
        return false;
    }


    // recursive coloring function
    Vector3d raytraceAPI::rayColor(Ray* ray, double t0, double t1, vector<Vector3d> lights, vector<Polygon*> polygons) {

        // if under five bounches & ray intersects with polygon
        if (ray->depth < DEPTH_CONST && findHit(ray, polygons)) {

            double lightRed = 0.0, lightGreen = 0.0, lightBlue = 0.0;
            double intensity = 1.0 / sqrt(lights.size());
            Vector3d returnVec(0, 0, 0);

            // for each light
            for (int r = 0; r < lights.size(); r++) {

                // building light ray object
                Vector3d lightNormal = (lights[r] - ray->intersections[0]->point);
                Ray* lightRay = new Ray(ray->intersections[0]->point + (EPSILON * lightNormal), lightNormal, EPSILON, ray->background);

                // if no shadow is found
                if (!findShadow(lightRay, polygons)) {

                    Polygon* polygon = ray->intersections[0]->object;
                    Vector3d candidate = polygon->color;


                    // building and normalizing vectors
                    Vector3d surfaceNormal = ray->intersections[0]->object->normal;
                    Vector3d eyeVector = -ray->direction;


                    surfaceNormal.normalize();
                    eyeVector.normalize();
                    lightNormal.normalize();

                    Vector3d halfway = (eyeVector + lightNormal);
                    halfway.normalize();

                    // creating refleciton direction via p87 in textbook
                    Vector3d reflect = ray->direction - (2 * (ray->direction.dot(surfaceNormal)) * surfaceNormal);

                    // creating reflection ray by offsetting source by epsilon as detailed in textbook p86-p87
                    Ray* reflection = new Ray(ray->intersections[0]->point + (EPSILON * reflect), reflect, ray->depth + 1, ray->background);

                    // from project description
                    double diffuse = getLightCoeff(surfaceNormal, lightNormal);
                    double specular = pow(getLightCoeff(surfaceNormal, halfway), polygon->shine);

                    lightRed += ((polygon->kd * candidate(0) * diffuse) + (polygon->ks * specular)) * intensity;
                    lightGreen += ((polygon->kd * candidate(1) * diffuse) + (polygon->ks * specular)) * intensity;
                    lightBlue += ((polygon->kd * candidate(2) * diffuse) + (polygon->ks * specular)) * intensity;

                    // translating local color to total color
                    returnVec(0) = lightRed; returnVec(1) = lightGreen; returnVec(2) = lightBlue;
                    // recursive call
                    returnVec += polygon->ks * rayColor(reflection, EPSILON, inf, lights, polygons);
                    delete reflection;
                }
                delete lightRay;
            }
            return returnVec;
        }
        else if (ray->depth < DEPTH_CONST) {
            // if no intersection, calculate with background color
            return ray->background;
        }
        else {
            // return zero after 5 bounces
            Vector3d nullV(0.0, 0.0, 0.0);
            return nullV;
        }
    }

    int raytraceAPI::openNFF(string nffFileName,
        Vector3d& backgroundColor,
        vector<Vector3d>* lights,
        vector<Polygon*>* polygons,
        Camera ** camera
        ) {

        double red, green, blue;
        double xFrom, yFrom, zFrom, xAt, yAt, zAt, xUp, yUp, zUp, angle, edges;
        double resolutionX, resolutionY;
        double hither;
        double redF, greenF, blueF, kd, ks, shine, t, index_of_refraction;
        char option = 'z';
        FILE* fp;
        fp = fopen(nffFileName.c_str(), "r");
        cout << "Parsing..." << endl;
        char str[101];
        char* buf[101];
        fgets(str, 100, fp);

        string line;
        while (!feof(fp)) {
          
            option = str[0];
            

            if (option == 'b') {
                sscanf(str, "%c %lf %lf %lf", &buf, &red, &green, &blue);
                backgroundColor << red, green, blue;
            }
            else if (option == 'v') {
                fgets(str, 100, fp);
                sscanf(str, "%4s %lf %lf %lf", &buf, &xFrom, &yFrom, &zFrom);

                fgets(str, 100, fp);
                sscanf(str, "%2s %lf %lf %lf", &buf, &xAt, &yAt, &zAt);

                fgets(str, 100, fp);
                sscanf(str, "%2s %lf %lf %lf", &buf, &xUp, &yUp, &zUp);

                fgets(str, 100, fp);
                sscanf(str, "%5s %lf", &buf, &angle);

                fgets(str, 100, fp);
                sscanf(str, "%6s %lf", &buf, &hither);

                fgets(str, 100, fp);
                sscanf(str, "%10s %lf %lf", &buf, &resolutionX, &resolutionY);

                *camera = new Camera(xFrom, yFrom, zFrom, xAt, yAt, zAt, xUp, yUp, zUp, angle, hither, resolutionX, resolutionY);
            }
            else if (option == 'f') {
                // TODO: seperate class for properties
                sscanf(str, "%c %lf %lf %lf %lf %lf %lf %lf %lf", &buf, &redF, &greenF, &blueF, &kd, &ks,
                    &shine, &t, &index_of_refraction);
            }
            else if (option == 'p') {
                sscanf(str, "%2s %lf", &buf, &edges);
                vector <Vector3d> vertices;
                Vector3d color(redF, greenF, blueF);
                for (int i = 0; i < edges; i++) {
                    Vector3d coordinates;
                    fgets(str, 100, fp);
                    sscanf(str, "%f %f %f", &coordinates(0), &coordinates(1), &coordinates(2));
                    vertices.push_back(coordinates);
                }
                polygons->push_back(new Polygon(edges, vertices, color, kd, ks, shine, index_of_refraction, t));
            }
            else if (option == 's') {
                // TODO: spheres
                double x, y, z, radius;
                sscanf(str, "%c %lf %lf %lf %lf", &buf, &x, &y, &z, &radius);
            }
            else if (option == 'l') {
                double x, y, z;
                sscanf(str, "%c %lf %lf %lf", &buf, &x, &y, &z);
                Vector3d light(x, y, z);
                lights->push_back(light);
            }
            else {
                cout << "Line unread " << endl;
            }
            fgets(str, 100, fp);
        }
        fclose(fp);
        return 0;
    }

    Ray* raytraceAPI::FireRayAtPixel(Camera* camera, double pixel_x, double pixel_y, Vector3d bgColor, vector<Vector3d> lights, vector<Vector3d> polygons) {
        cout << "Firing rays..." << endl;

        Vector3d pixel_s;
        camera->getPixelCoordinatesForPoint(pixel_x, pixel_x, pixel_s);


        Vector3d direction;
        direction = pixel_s[0] * camera->u + pixel_s[1] * camera->v + pixel_s[2] * camera->w;


        return new Ray(camera->from, direction, 0, bgColor);
    }

    vector<Vector3d> raytraceAPI::colorPixel(Ray* ray, vector<Vector3d> lights, vector<Polygon*> polygons) {

        double lightRed = 0.0, lightBlue = 0.0, lightGreen = 0.0;

        Vector3d pixelColor = rayColor(ray, 0, inf, lights, polygons);
        lightRed = pixelColor(0);
        lightGreen = pixelColor(1);
        lightBlue = pixelColor(2);

        vector<Vector3d> pixel;
        for (int x = 0; x < 3; x++) {
            // capping the light const at 1 (via prof bargteil on piazza)
            Vector3d color;
            color << min(lightRed, 1.0) * 255.0, min(lightGreen, 1.0) * 255.0, min(lightBlue, 1.0) * 255.0;
            pixel.push_back(color);
        }
        return pixel;
    }

    void raytraceAPI::deleteAllocated(vector<Polygon*> polygons, vector<Sphere*> spheres) {
        for (int i = 0; i < polygons.size(); i++) {
            delete polygons[i];
        }
        for (int i = 0; i < spheres.size(); i++) {
            delete spheres[i];
        }
        return;
    }


    int raytraceAPI::printImageFile(vector<vector<Vector3d > > pixels, Camera* camera) {
        computeTeapotFile(pixels, camera->resX, camera->resY);
        cout << "Completed Teapot File. View file 'teapot-seaman.ppm'  " << endl;
        return 0;
    }
