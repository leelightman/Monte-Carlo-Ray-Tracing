#ifndef CAMERA_H
#define CAMERA_H

#include <Eigen/Core>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include "Color.h"
using namespace std;
using namespace glm;
using namespace Eigen;

class Camera
{
public:
	vec3 eye;
	vec3 center;
	vec3 up;
	float fov;
	float aperture;
	float f;
	const int W;
	const int H;
	//mat4 V_P_inverse;

	Camera(const int w, const int h, const vec3 eye, const vec3 center, const vec3 up, const float fov, const float aperture, const float f);
	Ray initRay(const int i, const int j); // pixel loc
		
};

#endif