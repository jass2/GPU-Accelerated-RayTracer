/* 
*  Driver for teapot raytracing using OMP
*/


#include <iostream>
#include "Eigen/Dense"
#include "polygon.h"
#include "sphere.h"
#include "intersection.h"
#include "ray.h"
#include "raytraceAPI.h"
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include "Camera.h"
#include <omp.h>

using namespace std;

int main() {

	fflush(stdout);
	
	// make raytraceAPI class instance
	raytraceAPI rtapi;

	Camera* camera = NULL;
	Vector3d backgroundColor;
	vector < Vector3d > lights;
	vector < Polygon* > polygons;

	string file = "teapot-3.nff";

	// read in file 
	rtapi.openNFF(file, backgroundColor, &lights, &polygons, &camera);
	
	int resX = camera->resX;
	int resY = camera->resY;

	cout << "image " << file << " of size (" << resX << ", " << resY << ")" << endl;

	// output image pixels
	vector < vector < Vector3d > > pixels;

	// init the pixels cause vectors
	for (int i = 0; i < resX; i++) {
		vector<Vector3d> row;
		for (int j = 0; j < resY; j++) {
			Vector3d item;
			row.push_back(item);
		}
		pixels.push_back(row);
	}

	cout << "Firing Rays!" << endl;

	double start = omp_get_wtime();

	int i = 0;
	int j = 0;

	// for every pixel in the output image
//#pragma omp parallel for collapse(2) private(i, j) shared(pixels, backgroundColor, polygons, lights, camera) // action
	for (i = 0; i<resX; i++) {
		for (j = 0; j<resY; j++) {
			// shoot a ray to figure out where it ends up and get the color for it
			
			//Ray* FireRayAtPixel(Camera* camera, double pixel_x, double pixel_y, Vector3d bgColor, vector<Vector3d> lights);
			Ray* ray = rtapi.FireRayAtPixel(camera, i, j, backgroundColor, lights);
			Vector3d color = rtapi.colorPixel(ray, lights, polygons);

			// put color in right location in image
			pixels[i][j] = color;
		}
	}
	double end = omp_get_wtime();
	printf("Raytracing took %f seconds\n", end - start);

	// output file
	rtapi.printImageFile(pixels, camera);
	
	return 0;
}