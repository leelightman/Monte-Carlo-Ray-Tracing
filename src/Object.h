#ifndef OBJECT_H
#define OBJECT_H
#include <vector>
#include <Eigen/Dense>
#include <glm/glm.hpp>
#include "Color.h"
using namespace std;
using namespace glm;
using namespace Eigen;

class Object
{
private:
	const Material* material_;
public:
	Object(Material* material);
	virtual ~Object(){};

	virtual bool 	intersect(InterPt* id, Ray r) const = 0;
	Material 		material() const;
};


class Sphere : public Object
{
private:
	const vec3 POSITION_;
	const float RADIUS_;
public:
	Sphere(vec3 position, float radius, Material* material);

	bool intersect(InterPt* id, Ray r) const;
	vec3 get_pt_on_surface(float u, float v) const;
};

// P0, P1, and P2 defines a paralellogram which is the plane
class Plane : public Object
{
private:
	const vec3 P0_, P1_, P2_, V1_, V2_, NORMAL_;
	const float AREA_;
public:
	Plane(vec3 p0, vec3 p1, vec3 p2, Material* material);

	bool 	intersect(InterPt* id, Ray r) const;
	vec3 	get_pt_on_surface(float u, float v) const;
	float 	get_area() const;
	vec3 	get_normal() const;
	vec3 	get_first_tangent() const;
};

class Light
{
private:
	const Plane emitter_;
public:
	Light(
		vec3 p0,
		vec3 p1,
		vec3 p2,
		float flux, // Gets multiplied with color for total flux [Watts]
		Color color);
	
	bool 	intersect(LightInterPt* light_id, Ray r);
	vec3 	get_pt_on_surface(float u, float v);
	float 	get_area() const;
	vec3 	get_normal() const;
	Ray shoot_light_ray();
	const Color radiosity;
};

class Mesh : public Object
{
private:
	vector<vec3> positions_;
	vector<unsigned int> indices_;

	mat4 transform_;

public:
	Mesh(mat4 transform, string file_path, Material * material);
	
	bool intersect(InterPt* id, Ray r) const;
	vec3 get_min_pos() const;
	vec3 get_max_pos() const;
	mat4 get_transform() const;
	int	get_tri_num() const;
};

#endif