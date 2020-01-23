// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <omp.h>
#include <Eigen/Dense>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include <stdint.h>
#include <limits>
#include <cmath>
#include <ctime>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Camera.h"
#include "Scene.h"

// Shortcut to avoid Eigen:: and std::everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;
using namespace glm;
// Declare structure and classes


const double EPSILON = 0.0000001;
const double PHONG_COEFFICIENT = 50.0; 


void map_color_to_mat(Color specular, Color diffuse, MatrixXd &C_r, MatrixXd &C_g, MatrixXd &C_b,int i, int j){
    float gamma = 1/2.2;
    C_r(i,j) = clamp(pow(specular[0],gamma), 0.0f, 1.0f)+clamp(pow(diffuse[0],gamma), 0.0f, 1.0f);
    C_g(i,j) = clamp(pow(specular[1],gamma), 0.0f, 1.0f)+clamp(pow(diffuse[1],gamma), 0.0f, 1.0f);
    C_b(i,j) = clamp(pow(specular[2],gamma), 0.0f, 1.0f)+clamp(pow(diffuse[2],gamma), 0.0f, 1.0f);
}

Color ray_tracing(Color c, vec3 ray_dir, vec3 cam_plane_n){
    Color r;
    r  = c * dot(ray_dir, cam_plane_n);
    return r;
}
void mc_rt()
{
    cout << "Monte Carlo Ray Tracing. Please input your file name: ";
    string input_name;
    cin>> input_name;

    time_t time_start, time_end;
	time(&time_start);

	//const int W = 1024;
	//const int H = 768;
    //const int H = 1152;
	const int W = 100;
	const int H = 100;


	float progress = 0;
	cout << "Progress: " << endl;
	cout << progress << "\%" << endl;

    MatrixXd C_r = MatrixXd::Zero(W,H); // Red
    MatrixXd C_g = MatrixXd::Zero(W,H); // Green
    MatrixXd C_b = MatrixXd::Zero(W,H); // Blue
    MatrixXd A = MatrixXd::Ones(W,H); // Alpha

    const int SAMPLES = 100;
	
    Scene s;
    Camera c(W,H,vec3(0, 0, 2),vec3(0, 0, 0), vec3(0, 1, 0), M_PI / 2, 0.08,0.5);
	vec3 c_normal = normalize(c.center - c.eye);
	#pragma omp parallel for private(i)
    for (int i=0;i< C_r.rows(); i++)
	{
		// Parallellize the for loop with openMP.
		#pragma omp parallel for private(j)
		for (int j = 0; j< C_r.cols(); j++)
		{
            Color specular;
            Color diffuse;
            #pragma omp parallel for private(k)
			for (int k = 0; k < SAMPLES; k++)
			{
				//Ray r = c.initRay(i,(H - j - 1));
                Ray r = c.initRay(i,j); 
                specular += ray_tracing(s.trace_ray(r, 1), r.direction, c_normal);
                diffuse += ray_tracing(s.trace_ray(r, 2), r.direction, c_normal);
            }
            specular = specular / SAMPLES * (2 * M_PI);
            diffuse = diffuse/ SAMPLES * (2 * M_PI);
            //specular = specular / SAMPLES;
            //diffuse = diffuse/ SAMPLES;
			
            map_color_to_mat(specular,diffuse, C_r,C_g,C_b,i,j);
		}

		// Progress
        if(i%50 == 0){
            progress = (i+1) * 100 / float(W);
            cout << progress << "\%" << endl;
        }
	}

    if(1==1){
        time(&time_end);
        double ren_duration = difftime(time_end, time_start);
        int hr = ren_duration / (60 * 60);
        int min = (int(ren_duration) % (60 * 60)) / 60;
        int sec = int(ren_duration) % 60;
        string time_str =to_string(hr) + "h:"+ to_string(min) + "m:"+ to_string(sec) + "s";
        cout << "Total time: " << time_str << endl;
    }
    write_matrix_to_png(C_r,C_g,C_b,A,input_name+".png");
}

int main() {
    mc_rt();
	return 0;
}
