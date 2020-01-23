#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <map>
#include <random>
#include <Eigen/Dense>
#include <glm/glm.hpp>

#include "Color.h"
#include "Object.h"

using namespace std;
using namespace glm;
using namespace Eigen;

class Scene
{
private:
	random_device rd_;
	mt19937* gen_;
	uniform_real_distribution<float>* dis_;

	vector<Object*> objects_;
	vector<Light*> lights_;
	map<string, Material*> materials_;
	
  	// Normal path tracing for diffuse ray
	Color calc_diffuse(Ray r,int render_mode,InterPt id,int iteration);
	Color direct_diffuse(Ray r,int render_mode,InterPt id);
	Color global_diffuse(Ray r,int render_mode,InterPt id,int iteration);
	
	// Specular and refractive tracing
	Color trace_mirror_reflection(Ray r,int render_mode,InterPt id,int iteration);
	Color trace_refracted_ray(Ray r,int render_mode,InterPt id,int iteration,vec3 offset,bool inside);

	bool intersect(InterPt* id, Ray r);
	bool intersect_light(LightInterPt* light_id, Ray r);
public:
	Scene();
	enum RenderMode{BUFFER, SPECULAR, MONTE_CARLO};
	Color trace_ray(Ray r, int render_mode, int iteration = 0);

};

#endif // SCENE_H