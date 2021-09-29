#include <iostream>
#include <vector>
#include <Eigen>
#include <CL/cl.hpp>// main OpenCL include file 
#include "../../../483/Another/Project2/raytraceAPI.h"
#include "../../../483/Another/Project2/Camera.h"
#include "../polygon.h"
#include "../ray.h"
#include "../intersection.h"
#include <filesystem>
#include <chrono>



using namespace cl;
using namespace std;
namespace stdfs = filesystem;


cl_double4* cpu_output;
CommandQueue queue;
Kernel kernel;
Context context;
Program program;
Buffer cl_output, cl_polygons, cl_lights, cl_camera;
Device device;

struct st_Camera {
	cl_double angle;
	cl_double hither;
	cl_double resX;
	cl_double resY;
	cl_double top;
	cl_double right;
	cl_double bottom;
	cl_double left;
	cl_double3 from;
	cl_double3 at;
	cl_double3 up;
	cl_double3 u;
	cl_double3 v;
	cl_double3 w;
	cl_double3 background;
};

struct st_Polygon {
	cl_double numEdges;
	cl_double kd;
	cl_double ks;
	cl_double shine;
	cl_double refraction;
	cl_double kt;
	cl_double d1;
	cl_double d2;
	cl_double3 v1;
	cl_double3 v2;
	cl_double3 v3;
	cl_double3 v4;
	cl_double3 color;
	cl_double3 normal;
};

void pickPlatform(Platform& platform, const vector<Platform>& platforms) {

	if (platforms.size() == 1) platform = platforms[0];
	else {
		int input = 0;
		cout << "\nChoose an OpenCL platform: ";
		cin >> input;

		// handle incorrect user input
		while (input < 1 || input > platforms.size()) {
			cin.clear(); //clear errors/bad flags on cin
			cin.ignore(cin.rdbuf()->in_avail(), '\n'); // ignores exact number of chars in cin buffer
			cout << "No such option. Choose an OpenCL platform: ";
			cin >> input;
		}
		platform = platforms[input - 1];
	}
}

void pickDevice(Device& device, const vector<Device>& devices) {

	if (devices.size() == 1) device = devices[0];
	else {
		int input = 0;
		cout << "\nChoose an OpenCL device: ";
		cin >> input;

		// handle incorrect user input
		while (input < 1 || input > devices.size()) {
			cin.clear(); //clear errors/bad flags on cin
			cin.ignore(cin.rdbuf()->in_avail(), '\n'); // ignores exact number of chars in cin buffer
			cout << "No such option. Choose an OpenCL device: ";
			cin >> input;
		}
		device = devices[input - 1];
	}
}

void printErrorLog(const Program& program, const Device& device) {

	// Get the error log and print to console
	string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
	cerr << "Build log:" << std::endl << buildlog << std::endl;

	// Print the error log to a file
	FILE* log = fopen("errorlog.txt", "w");
	fprintf(log, "%s\n", buildlog);
	cout << "Error log saved in 'errorlog.txt'" << endl;
	system("PAUSE");
	exit(1);
}

void initOpenCL()
{
	// Get all available OpenCL platforms (e.g. AMD OpenCL, Nvidia CUDA, Intel OpenCL)
	vector<Platform> platforms;
	Platform::get(&platforms);
	cout << "Available OpenCL platforms : " << endl << endl;
	for (int i = 0; i < platforms.size(); i++)
		cout << "\t" << i + 1 << ": " << platforms[i].getInfo<CL_PLATFORM_NAME>() << endl;

	// Pick one platform
	Platform platform;
	pickPlatform(platform, platforms);
	cout << "\nUsing OpenCL platform: \t" << platform.getInfo<CL_PLATFORM_NAME>() << endl;

	// Get available OpenCL devices on platform
	vector<Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

	cout << "Available OpenCL devices on this platform: " << endl << endl;
	for (int i = 0; i < devices.size(); i++) {
		cout << "\t" << i + 1 << ": " << devices[i].getInfo<CL_DEVICE_NAME>() << endl;
		cout << "\t\tMax compute units: " << devices[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << endl;
		cout << "\t\tMax work group size: " << devices[i].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << endl << endl;
	}

	// Pick one device
	pickDevice(device, devices);
	cout << "\nUsing OpenCL device: \t" << device.getInfo<CL_DEVICE_NAME>() << endl;
	cout << "\t\t\tMax compute units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << endl;
	cout << "\t\t\tMax work group size: " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << endl;

	// Create an OpenCL context and command queue on that device.
	context = Context(device);
	queue = CommandQueue(context, device);

	// Convert the OpenCL source code to a string
	string source;
	ifstream file("kernel.cl", ios::binary);
	if (!file) {
		cout << "\nNo OpenCL file found!" << endl << "Exiting..." << endl;
		system("PAUSE");
		exit(1);
	}
	while (!file.eof()) {
		char line[256];
		file.getline(line, 255);
		source += line;
	}

	const char* kernel_source = source.c_str();

	// Create an OpenCL program by performing runtime source compilation for the chosen device
	program = Program(context, kernel_source);
	cl_int result = program.build({ device });
	if (result) cout << "Error during compilation OpenCL code!!!\n (" << result << ")" << endl;
	if (result == CL_BUILD_PROGRAM_FAILURE) printErrorLog(program, device);

	// Create a kernel (entry point in the OpenCL source program)
	kernel = Kernel(program, "render_kernel");
}

void cleanUp() {
	delete cpu_output;
}

inline double clamp(double x) { return x < 0.0f ? 0.0f : x > 1.0f ? 1.0f : x; }

// convert RGB float in range [0,1] to int in range [0, 255] and perform gamma correction
inline int toInt(double x) { return int(clamp(x) * 255 + .5); }

void saveImage() {
	// write image to PPM file, a very simple image file format
	// PPM files can be opened with IrfanView (download at www.irfanview.com) or GIMP
	FILE* f = fopen("opencl_raytracer.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n", 512, 512, 255);

	// loop over all pixels, write RGB values
	for (int i = 0; i < 512; i++) {
		for (int j = 0; j < 512; j++) {
			fprintf(f, "%d %d %d ",
				toInt(cpu_output[i * 512 + 511 - j].s[0]),
				toInt(cpu_output[i * 512 + 511 - j].s[1]),
				toInt(cpu_output[i * 512 + 511 - j].s[2]));
		}
	}
}

cl_double3 convertVectTodouble(Vector3d bg) {
	cl_double3 color = { (double)bg[0], (double)bg[1], (double)bg[2] };
	return color;
}

int main()
{

	raytraceAPI rtapi;
	Vector3d backgroundColor;
	vector<Vector3d> lights;
	vector<Polygon*> polygons;
	vector<vector<Vector3d> > pixels;
	Camera * camera = NULL;
	
	struct st_Camera st_cam;

	string filename;
	stdfs::path cwd;
	cwd = stdfs::current_path();
	cout << cwd.string() << endl;
	cout << endl << "Select an input filename: " << endl;
	cin >> filename;
	rtapi.openNFF(filename, backgroundColor, &lights, &polygons, &camera);
	
	int polys = (int) polygons.size();
	int numItems = (int)(camera->resX * camera->resY);
	int numLights = (int)lights.size();
	cl_double3 bc = convertVectTodouble(backgroundColor);

	st_cam.angle = (double)camera->angle;
	st_cam.hither = (double)camera->hither;
	st_cam.resX = (double)camera->resX;
	st_cam.resY = (double)camera->resY;
	st_cam.top = (double)camera->top;
	st_cam.right = (double)camera->right;
	st_cam.bottom = (double)camera->bottom;
	st_cam.left = (double)camera->left;
	st_cam.from = { (double)camera->from[0], (double)camera->from[1], (double)camera->from[2] };
	st_cam.at = { (double)camera->at[0], (double)camera->at[1], (double)camera->at[2] };
	st_cam.up = { (double)camera->up[0], (double)camera->up[1], (double)camera->up[2] };
	st_cam.u = { (double)camera->u[0], (double)camera->u[1], (double)camera->u[2] };
	st_cam.v = { (double)camera->v[0], (double)camera->v[1], (double)camera->v[2] };
	st_cam.w = { (double)camera->w[0], (double)camera->w[1], (double)camera->w[2] };
	st_cam.background = bc;

	cpu_output = new cl_double3[numItems];

	cl_double3 st_lights[2];
	st_Polygon st_polygons[600];

	for (int i = 0; i < polys; i++) {
		st_polygons[i].numEdges = (double)polygons[i]->numEdges;
		st_polygons[i].v1 = { (double)polygons[i]->vertices[0][0], (double)polygons[i]->vertices[0][1], (double)polygons[i]->vertices[0][2] };
		st_polygons[i].v2 = { (double)polygons[i]->vertices[1][0], (double)polygons[i]->vertices[1][1], (double)polygons[i]->vertices[1][2] };
		st_polygons[i].v3 = { (double)polygons[i]->vertices[2][0], (double)polygons[i]->vertices[2][1], (double)polygons[i]->vertices[2][2] };
		if (polygons[i]->numEdges == 4)
			st_polygons[i].v4 = { (double)polygons[i]->vertices[3][0], (double)polygons[i]->vertices[3][1], (double)polygons[i]->vertices[3][2] };
		else
			st_polygons[i].v4 = { (double)0.0f, (double)0.0f, (double)0.0f };
		st_polygons[i].kd = (double)polygons[i]->kd;
		st_polygons[i].ks = (double)polygons[i]->ks;
		st_polygons[i].shine = (double)polygons[i]->shine;
		st_polygons[i].refraction = (double)polygons[i]->refraction;
		st_polygons[i].kt = (double)polygons[i]->kt;
		st_polygons[i].normal = { (double)polygons[i]->normal[0], (double)polygons[i]->normal[1], (double)polygons[i]->normal[2] };
		st_polygons[i].color = { (double)polygons[i]->color[0], (double)polygons[i]->color[1], (double)polygons[i]->color[2] };
	}

	for (int i = 0; i < numLights; i++) {
		st_lights[i] = { lights[i][0], lights[i][1], lights[i][2] };
	}

	// initialise OpenCL
	initOpenCL();

	auto start = chrono::high_resolution_clock::now();

	// Create image buffer on the OpenCL device
	cl_output = Buffer(context, CL_MEM_WRITE_ONLY, numItems * sizeof(cl_double3));
	cl_polygons = Buffer(context, CL_MEM_READ_ONLY, polys * sizeof(st_Polygon));
	cl_lights = Buffer(context, CL_MEM_READ_ONLY, numLights * sizeof(cl_double3));
	cl_camera = Buffer(context, CL_MEM_READ_ONLY, sizeof(st_Camera));

	queue.enqueueWriteBuffer(cl_polygons, CL_TRUE, 0, (int)polygons.size() * sizeof(st_Polygon), st_polygons);
	queue.enqueueWriteBuffer(cl_lights, CL_TRUE, 0, (int)lights.size() * sizeof(cl_double3), st_lights);
	queue.enqueueWriteBuffer(cl_camera, CL_TRUE, 0, sizeof(st_Camera), &st_cam);

 	kernel.setArg(0, cl_output);
	kernel.setArg(1, (int)camera->resY);
	kernel.setArg(2, (int)camera->resX);
	kernel.setArg(3, cl_camera);
	kernel.setArg(4, cl_polygons);
	kernel.setArg(5, cl_lights);
	kernel.setArg(6, numLights);
	kernel.setArg(7, polys);

	// every pixel in the image has its own thread or "work item",
	// so the total amount of work items equals the number of pixels
	std::size_t global_work_size = (512 * 512);
	std::size_t local_work_size = kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
	cout << "Kernel work group size: " << local_work_size << endl;

	// Ensure the global work size is a multiple of local work size
	if (global_work_size % local_work_size != 0)
		global_work_size = (global_work_size / local_work_size + 1) * local_work_size;

	// launch the kernel
	queue.enqueueNDRangeKernel(kernel, NULL, global_work_size, local_work_size);
	queue.finish();

	// read and copy OpenCL output to CPU
	queue.enqueueReadBuffer(cl_output, CL_TRUE, 0, 512 * 512 * sizeof(cl_double3), cpu_output);

	// save image to PPM format
	saveImage();

	cout << "Rendering done!\nSaved image to 'teapot.ppm'" << endl;
	// record and display runtime
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
	cout << duration.count() << "microseconds " << endl;

	// release memory
	cleanUp();

	return 0;
}