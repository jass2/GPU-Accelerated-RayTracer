//
// Created by Jason Seaman on 2021-05-02.
//

#ifndef JASS2_CAMERA_H
#define JASS2_CAMERA_H

#include <Eigen>
#include <cmath>
#include <vector>

using namespace std;

typedef Eigen::Matrix <float, 3, 1 > Vector3d;


// proj1 funciton


class Camera {

public:
    Vector3d from;
    Vector3d at;
    Vector3d up;

    Vector3d w;
    Vector3d u;
    Vector3d v;

    double angle;
    double hither;
    double resX;
    double resY;

    double top;
    double right;
    double bottom;
    double left;

    Camera() {};

    Camera(double xFrom, double yFrom, double zFrom, double xAt, double yAt, double zAt, double xUp, double yUp, double zUp, double viewAngle, double _hither, double _resolutionX, double _resolutionY) {
        
        from << xFrom, yFrom, zFrom;
        at << xAt, yAt, zAt;
        up << xUp, yUp, zUp;
        
        angle = getRadians(viewAngle);

        top = tan(angle / 2);
        right = top;
        bottom = -top;
        left = -right;

        w = from - at;
        w.normalize();
       
        u = up.cross(w);
        u.normalize();
        
        // implicitly normalized
        v = w.cross(u);

        resX = _resolutionX;
        resY = _resolutionY;
        hither = _hither;
    }

    double getRadians(double value) {
        return value * (3.14159262 / 180);
    }

    void getPixelCoordinatesForPoint(double x, double y, Vector3d& p_s) {

        double u_s = left + (right - left) * (x + 0.5) / resX;
        double v_s = bottom + (top - bottom) * (y + 0.5) / resY;
        double w_s = -1;

        p_s << u_s, v_s, w_s;
        return;
    }

    ~Camera() {}

};




#endif
