/* OpenCL ray tracing tutorial by Sam Lapere, 2016
http://raytracey.blogspot.com */

typedef struct st_Polygon {
	double numEdges;
	double kd;
	double ks;
	double shine;
	double refraction;
	double kt;
	double d1;
	double d2;
	double3 v1;
	double3 v2;
	double3 v3;
	double3 v4;
	double3 color;
	double3 normal;
} st_Polygon;

typedef struct Ray {
	double3 origin;
	double3 dir;
	double3 eye;
	double depth;
	double3 background;
	double3 closest;
	double distance;
	double3 color;
	st_Polygon intersectObject;
} Ray;

typedef struct st_Camera {
	double angle;
	double hither;
	double resX;
	double resY;
	double top;
	double right;
	double bottom;
	double left;
	double3 from;
	double3 at;
	double3 up;
	double3 u;
	double3 v;
	double3 w;
	double3 background;
} st_Camera;

__constant static double inf = 1e5f;
__constant static int DEPTH_CONST = 3.0;
__constant static double EPSILON = 0.001f;

void printVec(double3 item) {
	printf("vert = %lf, %lf, %lf\n", item.x, item.y, item.z);
}

double3 getPixelCoordinatesForPoint(__constant st_Camera * camera, const double x, const double y) {

	double u_s = camera->left + (camera->right - camera->left) * (x + 0.5f) / camera->resX;
	double v_s = camera->bottom + (camera->top - camera->bottom) * (y + 0.5f) / camera->resY;
	double w_s = -1.0f;
	double3 pixel = { u_s, v_s, w_s };
	return pixel;
}

struct Ray FireRayAtPixel(__constant st_Camera * camera, const double pixel_x, const double pixel_y) {
	
	double3 pixel_s = getPixelCoordinatesForPoint(camera, pixel_x, pixel_y);
	Ray ray;
	ray.origin = camera->from;
	ray.dir = pixel_s.z * camera->w + pixel_s.x * camera->u + pixel_s.y * camera->v;
	ray.depth = 0;
	ray.background = camera->background;
	ray.eye = camera->from;
	ray.color = camera->background;
	ray.closest = (double3){ 0.0f, 0.0f, 0.0f };
	return ray;
}

bool intersect_poly(__constant st_Polygon* polygon, int four, const Ray* ray, const double t0, double t1, double3* pointHit) {
	double M, t, gamma, beta;
	double3 abc, def, ghi, jkl;
	if (!four) {
		abc = polygon->v1 - polygon->v2;
		def = polygon->v1 - polygon->v3;
	}
	else {
		abc = polygon->v1 - polygon->v3;
		def = polygon->v1 - polygon->v4;
	}
	ghi = ray->dir;
	jkl = (polygon->v1 - ray->origin);

		M = abc.x * (def.y * ghi.z - ghi.y * def.z) +
			abc.y * (def.z * ghi.x - ghi.z * def.x) +
			abc.z * (def.x * ghi.y - ghi.x * def.y);

		t = -1.0f * (
			(def.z * (abc.x * jkl.y - jkl.x * abc.y)) +
			(def.y * (abc.z * jkl.x - jkl.z * abc.x)) +
			(def.x * (abc.y * jkl.z - jkl.y * abc.z))
			) / M;

		if (t < t0 || t > t1) {
			return false;
		}

		gamma = (
			(ghi.z * (abc.x * jkl.y - jkl.x * abc.y)) +
			(ghi.y * (abc.z * jkl.x - jkl.z * abc.x)) +
			(ghi.x * (abc.y * jkl.z - jkl.y * abc.z))
			) / M;

		if (gamma < 0.0f || gamma > 1.0f) {
			return false;
		}

		beta = (
			(jkl.x * (def.y * ghi.z - def.z * ghi.y)) +
			(jkl.y * (def.z * ghi.x - def.x * ghi.z)) +
			(jkl.z * (def.x * ghi.y - def.y * ghi.x))
			) / M;

		if (beta < 0.0f || (beta > 1.0f - gamma)) {
			return false;
		}
	
	if (!four) {
		*pointHit = polygon->v1 + beta * (polygon->v2 - polygon->v1) + gamma * (polygon->v3 - polygon->v1);
	}
	else {
		*pointHit = polygon->v1 + beta * (polygon->v3 - polygon->v1) + gamma * (polygon->v4 - polygon->v1);
	}
	return true;
}

// test for first intersection between point of intersection and light source
bool findShadow(Ray * light, __constant st_Polygon * polygons, const int numPolygons) {
	double3 pointShadow;
	bool hit = false;
	double closestDist = inf;
	for (int i = 0; i < numPolygons; i++) {
		if (intersect_poly(&polygons[i], 0, light, 0.0f, inf, &pointShadow)) {
			hit = true;
			i = numPolygons;
		}
		else if (polygons[i].numEdges > 3.001f){
			if (intersect_poly(&polygons[i], 1, light, 0.0f, inf, &pointShadow)) {
				hit = true;
				i = numPolygons;
			}
		}
		else {
			continue;
		}
	}
	return hit;
}


bool testSubPolys(Ray * ray, __constant st_Polygon * polygons, const int numPolygons, const double t0, const double t1) {
	bool intersect = false;
	double closestDist = t1;
	double3 closestIntersection;
	double3 closestColor;
	st_Polygon intersectObject;
	double3 pointOfIntersect;
	double distFromIntersect;
	for (int i = 0; i < numPolygons; i++) {
		if (intersect_poly(&polygons[i], 0, ray, t0, closestDist, &pointOfIntersect)) {
			distFromIntersect = distance(pointOfIntersect, ray->eye);
			if (min(distFromIntersect, closestDist) < closestDist) {
				closestDist = distFromIntersect;
				closestIntersection = pointOfIntersect;
				intersectObject = polygons[i];
				closestColor = polygons[i].color;
			}
			intersect = true;
		}
		if (polygons[i].numEdges > 3.001f){
			if (intersect_poly(&polygons[i], 1, ray, t0, closestDist, &pointOfIntersect)) {
				distFromIntersect = distance(pointOfIntersect, ray->eye);
				if (min(distFromIntersect, closestDist) < closestDist) {
					closestDist = distFromIntersect;
					closestIntersection = pointOfIntersect;
					intersectObject = polygons[i];
					closestColor = polygons[i].color;
				}
				intersect = true;
			}
		}
	}
	if (intersect) {
		ray->intersectObject = intersectObject;
		ray->closest = closestIntersection;
		ray->distance = closestDist;		
		ray->color = closestColor;
	}
	return intersect;
}


double getLightCoeff(const double3 normal, const double3 dottedWith) {
	double num = dot(normal, dottedWith);
	if (num > 0.0f) {
		return num;
	}
	else {
		return 0.0f;
	}
}

double3 colorRay(Ray * ray,const double t0, const double t1, __constant double3* lights, __constant st_Polygon * polygons, const int numLights, const int numPolygons) {

	double intensity = 1.0f / sqrt((double)numLights);
	double3 lightColor = { 0.0f, 0.0f, 0.0f };
	double lastKs = 1.0f;
	Ray lightRay;
	for (int i = 0; i < DEPTH_CONST; i++) {
		if (!testSubPolys(ray, polygons, numPolygons, t0, t1)) {
			return lightColor += lastKs * ray->background;
		}
		else {
			st_Polygon polygon = ray->intersectObject;
			/* compute the surface normal and flip it if necessary to face the incoming ray */
			double3 surfaceNormal = normalize(dot(polygon.normal, ray->eye) > 0.0f ? polygon.normal : polygon.normal * (-1.0f));
			for (int r = 0; r < numLights; r++) {
				lightRay.dir = normalize(lights[r] - ray->closest);
				lightRay.origin = ray->closest + (EPSILON * lightRay.dir);
				lightRay.depth = EPSILON;

				if (!findShadow(&lightRay, polygons, numPolygons)) {
					double3 lightNormal = lightRay.dir;
					double3 eyeVector = -1.0f * ray->dir;
					double3 halfway = normalize(eyeVector + lightNormal);
					double diffuse = getLightCoeff(surfaceNormal, lightNormal);
					double specular = pow(getLightCoeff(surfaceNormal, halfway), polygon.shine);
					lightColor += lastKs * (polygon.kd * polygon.color * diffuse + polygon.ks * specular) * intensity;
				}
			}
			ray->dir = normalize(ray->dir - 2.0f * dot(ray->dir, surfaceNormal) * surfaceNormal);
			ray->origin = ray->closest + EPSILON * ray->dir;
			lastKs = polygon.ks;
		}
	}
	return lightColor;
}

__kernel void render_kernel(
	__global double3 * output,
	int width,
	int height, 
	__constant st_Camera * camera,
	__constant st_Polygon * polygons,
	__constant double3 * lights,
	int numLights,
	int numPolygons)
{

	const int work_item_id = get_global_id(0);		/* the unique global id of the work item for the current pixel */
	
	int x_coord = work_item_id % width;					/* x-coordinate of the pixel */
	int y_coord = work_item_id / width;					/* y-coordinate of the pixel */

	/*create a camera ray */
	Ray camray = FireRayAtPixel(camera, x_coord, y_coord);
	
	/*get the color of pixel*/
	double3 color = colorRay(&camray, 0.0f, inf, lights, polygons, numLights, numPolygons);
	
	/*write the color to the output buffer*/
	output[(width * height) - work_item_id - 1] = color;
}