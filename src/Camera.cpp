#include "Camera.h"
#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <math.h>

using namespace glm;
using namespace std;
using namespace Eigen;

Camera::Camera (int w, int h, vec3 eye, vec3 center, vec3 up, float fov, float aperture, float f) :  W(w), H(h), eye(eye), center(center), up(up), fov(fov), aperture(aperture), f(f)
{
	// View and perspective matrices are used in the unProject() function
	//mat4 V = lookAt(eye, center, up);
	float aspect = float(W) / H;
	//mat4 P = perspective(fov, aspect, 0.1f, 100.0f);
	//V_P_inverse = inverse(V * P);
}

Ray Camera::initRay(int i, int j)
{
	Ray r;
	if (i < 0 || i > W - 1 || j < 0 || j > H - 1)
	{
		cout << "error---" << endl;
		r.origin = vec3(0);
		r.direction = vec3(0);
	}
	else
	{
		// https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range
		// generate randon number in range [-0.5,0.5)
		random_device rand_dev;
		mt19937 generator(rand_dev());
		uniform_real_distribution<float> distr(-0.5, 0.5);
		uniform_real_distribution<float> distr2(-0.04, 0.04);

		// randomize displacement for anti-aliasing
		float x_rand= distr(generator);
		//cout<<"x_displacement: "<<x_displacement<<endl;
		float y_rand = distr(generator);
		
		//x_rand= 0;
		//y_rand = 0;
		/*
		// View and perspective matrices are used in the unProject() function
		mat4 V = lookAt(eye, center, up);
		float aspect = float(W) / H;
		mat4 P = perspective(fov, aspect, 0.1f, 100.0f);
		*/

		//vec4 ray_origin_4 = VP_inv * vec4(((i + x_rand) / W - 0.5) * 2, ((j + y_rand) / H - 0.5) * 2, 1, 1 ); 
		//vec4 ray_end_4 = VP_inv * vec4(((i + x_rand) / W - 0.5) * 2, ((j + y_rand) / H - 0.5) * 2, -1, 1 );
		//vec3 ray_origin = vec3(ray_origin_4) * ray_origin_4.w;
		//cout<<"ray_origin_4.w: "<<ray_origin_4.w<<endl;
		//vec3 ray_end = vec3(ray_end_4) * ray_end_4.w;
		//vec3 ray_dir = normalize(ray_end - ray_origin);
		vec3 origin(-1,1,1);
		vec3 ray_origin = eye;
		vec3 x_displacement(2.0/W,0,0);
     	vec3 y_displacement(0,-2.0/H,0); 
		vec3 ray_dir = normalize(origin + double(i)*x_displacement+double(j)*y_displacement-ray_origin);


		float x_dop = distr2(generator);
		float y_dop = sqrt((aperture/2)*(aperture/2) - x_dop*x_dop);
		//cout<<"ray_dir: "<<ray_dir[0]<<" "<<ray_dir[1]<<" "<<ray_dir[2]<<endl;
		vec3 new_focal_aim = eye + 2*ray_dir;
		//cout<<"new_focal_aim: "<<new_focal_aim[0]<< " "<<new_focal_aim[1]<< " "<<new_focal_aim[2]<<endl;
		//vec3 new_focal_aim (x_dop, y_dop, eye[2]-0.3);
		vec3 new_eye (eye[0]+x_dop, eye[1]+y_dop, eye[2]);
		vec3 new_dir = normalize(new_focal_aim-new_eye);
		//cout<<ray_dir[0]<<" "<<ray_dir[1]<<" "<<ray_dir[1]<<" "<<endl;
		//vec3 ray_dir = vec3(0,0,-1);

		//r.origin = eye;
		r.origin = new_eye;
		r.direction = new_dir;
		//r.direction = ray_dir;
		r.radiance = Color::white();
		r.material = Material::air();
	}
	return r;
}