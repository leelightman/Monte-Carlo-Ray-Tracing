#include "Scene.h"
#include "Color.h"
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <random>
using namespace std;
using namespace Eigen;
using namespace glm;
// --- Scene class functions --- //

Scene::Scene ()
{
	gen_ = new mt19937(rd_());
	dis_ = new uniform_real_distribution<float>(0, 1);

	Color spec;
	Color dif;

	Material* blue = new Material;
	dif[0] = 0.2;
	dif[1] = 0.2;
	dif[2] = 0.9;
	spec[0]=1;
	spec[1]=1;
	spec[2]=1;
    blue->color_diffuse = dif;
    blue->color_specular = spec;
    blue->reflectance = 1;
    blue->specular_reflectance = 0; 
    blue->transmissivity = 0;
    blue->refraction_index = 1;
	blue->name = "blue";

	Material* white = new Material;
	dif[0] = 1;
	dif[1] = 1;
    white->color_diffuse = dif;
    white->color_specular = spec;
    white->reflectance = 0.8;
    white->specular_reflectance = 0; 
    white->transmissivity = 0;
    white->refraction_index = 1;

	Material* gray = new Material;
	dif[0] = 0.5;
	dif[1] = 0.5;
	dif[2] = 0.5;
    gray->color_diffuse = dif;
    gray->color_specular = spec;
    gray->reflectance = 0.8;
    gray->specular_reflectance = 0; 
    gray->transmissivity = 0;
    gray->refraction_index = 1;

	Material* glass = new Material;
	dif[0] =1;
	dif[1] =1;
	dif[2] =1;
    glass->color_diffuse = dif;
    glass->color_specular = spec;
    glass->reflectance = 1;
    glass->specular_reflectance = 1; 
    glass->transmissivity = 1;
    glass->refraction_index = 1.6;
	
	Material* mirror = new Material;
	dif[0] =0.2;
	dif[1] =0.2;
	dif[2] =0.8;
	spec[1]=0.9;
	spec[2]=0.8;
    mirror->color_diffuse = dif;
    mirror->color_specular = spec;
    mirror->reflectance = 1;
    mirror->specular_reflectance = 1; 
    mirror->transmissivity = 0;
    mirror->refraction_index = 1;

	this->materials_.insert(pair<string, Material*>("blue", blue));
	this->materials_.insert(pair<string, Material*>("glass", glass));
	this->materials_.insert(pair<string, Material*>("gray", gray));
	this->materials_.insert(pair<string, Material*>("mirror", mirror));


	Object* back = new Plane(vec3(-1.5,-1,-1), vec3(1.5,-1,-1), vec3(-1.5,1.5,-1), this->materials_["gray"]);
	this->objects_.push_back(back);
	Object* left = new Plane(vec3(-1.5,-1,-1), vec3(-1.5,1.5,-1), vec3(-1.5,-1,6), this->materials_["gray"]);
	this->objects_.push_back(left);
	Object* right = new Plane(vec3(1.5,-1,-1), vec3(1.5,-1,6), vec3(1.5,1.5,-1), this->materials_["gray"]);
	this->objects_.push_back(right);
	Object* top = new Plane(vec3(-1.5,1.5,-1), vec3(1.5,1.5,-1), vec3(-1.5,1.5,6), this->materials_["gray"]);
	this->objects_.push_back(top);
	Object* bottom = new Plane(vec3(-1.5,-1,-1), vec3(-1.5,-1,6), vec3(1.5,-1,-1), this->materials_["gray"]);
	this->objects_.push_back(bottom);

	Object* ball1 = new Sphere(vec3(0.9,-0.7,-0.7), 0.3, this->materials_["glass"]);
	this->objects_.push_back(ball1);
	Object* ball2 = new Sphere(vec3(-0.9,-0.65,0.7), 0.3, this->materials_["mirror"]);
	this->objects_.push_back(ball2);
	Object* ball3 = new Sphere(vec3(0,-0.7,0.1), 0.35, this->materials_["blue"]);
	this->objects_.push_back(ball3);

	mat4 transform = scale(mat4(), vec3(0.4, 0.4, 0.4));
	Object* cube1 = new Mesh(transform,"../data/cube_tri.off",this->materials_["blue"]);
	//this->objects_.push_back(cube1);
	/*
	Object* cube_back = new Plane(vec3(-.5,-.5,0), vec3(.5,-.5,0), vec3(-.5,.5,0), this->materials_["gray"]);
	this->objects_.push_back(cube_back);
	Object* cube_left = new Plane(vec3(-.5,-.5,0), vec3(-.5,.5,0), vec3(-.5,-.5,1), this->materials_["gray"]);
	this->objects_.push_back(cube_left);
	Object* cube_right = new Plane(vec3(.5,-.5,0), vec3(.5,.5,0), vec3(.5,.5,1), this->materials_["gray"]);
	this->objects_.push_back(cube_right);
	Object* cube_top = new Plane(vec3(-.5,.5,0), vec3(.5,.5,0), vec3(.5,.5,1), this->materials_["gray"]);
	this->objects_.push_back(cube_top);
	Object* cube_bottom = new Plane(vec3(-.5,-.5,0), vec3(-.5,-.5,1), vec3(.5,-.5,0), this->materials_["gray"]);
	this->objects_.push_back(cube_bottom);
	Object* cube_front = new Plane(vec3(-.5,-.5,1), vec3(-.5,.5,1), vec3(.5,-.5,1), this->materials_["gray"]);
	this->objects_.push_back(cube_front);
	*/


	Color color;
	color[0] = 1;
	color[1] = 1;
	color[2] = 1;
	Light* light = new Light(vec3(-0.5,1.49999,-0.5), vec3(0.5,1.49999,-0.5), vec3(-0.5,1.49999,0.5), 12, color);
	this->lights_.push_back(light);
}

//recursive
bool Scene::intersect(InterPt* ip, Ray r)
{
	InterPt ip_min_t;
	ip_min_t.t = 1000; // initiate t

	Object* intersecting_obj = NULL;
	for (int i = 0; i < objects_.size(); ++i)
	{
		InterPt ip_local;
		if (objects_[i]->intersect(&ip_local,r) && ip_local.t < ip_min_t.t)
		{
			ip_min_t = ip_local;
			intersecting_obj = objects_[i];
		}
	}
	if (intersecting_obj)
	{
		*ip = ip_min_t;
		return true;
	}
	return false;
}

bool Scene::intersect_light(LightInterPt* light_ip, Ray r)
{
	LightInterPt light_ip_min_t;
	light_ip_min_t.t = 1000; 
	Light* intersecting_light = NULL;
	LightInterPt ip_local;
	if (lights_[0]->intersect(&ip_local,r) && ip_local.t < light_ip_min_t.t)
	{
		light_ip_min_t = ip_local;
		intersecting_light = lights_[0];
	}
	
	if (intersecting_light)
	{
		InterPt ip_min_t;
		ip_min_t.t = 1000;
		Object* intersecting_obj = NULL;
		for (int i = 0; i < objects_.size(); ++i)
		{
			InterPt ip_local;
			if (objects_[i]->intersect(&ip_local,r) && ip_local.t < ip_min_t.t)
			{
				ip_min_t = ip_local;
				intersecting_obj = objects_[i];
			}
		}
		if (intersecting_obj && ip_min_t.t < light_ip_min_t.t)
		{
			return false;
		}
		else
		{
			*light_ip = light_ip_min_t;
			return true;	
		}
	}
	return false;
}

// Diffuse = local illumination (shadow rays) + indirect
Color Scene::calc_diffuse(Ray r,int render_mode,InterPt ip,int iteration)
{
	r.has_intersected = true;
	Color dif;
	Color direct_dif = direct_diffuse(r, render_mode, ip);
	Color global_dif = global_diffuse(r, render_mode, ip, iteration);
	// global illumination with mc
	dif = direct_dif + global_dif;
	return dif;
}

Color Scene::direct_diffuse(Ray r,int render_mode,InterPt ip)
{
	Color direct;
	const double PHONG_COEFFICIENT = 50.0;
	// divide area light into n_samples
	static const int n_samples = 1;

	for (int j = 0; j < n_samples; ++j)
	{
		Ray shadow_ray = r;
		vec3 differance = lights_[0]->get_pt_on_surface((*dis_)(*gen_),(*dis_)(*gen_)) - shadow_ray.origin;
		shadow_ray.direction = normalize(differance);
		Color brdf; // Dependent on inclination and azimuth
		Color albedo = ip.material.color_diffuse * ip.material.reflectance;
		Color albedo_spec = ip.material.color_specular * ip.material.reflectance;
		float cos_theta = dot(shadow_ray.direction, ip.normal);
		LightInterPt shadow_ray_ip;

		//brdf = calc_lambertian_BRDF(-r.direction,shadow_ray.direction,ip.normal,albedo);
		brdf = calc_lambertian_BRDF(albedo);
		//if(ip.material.name == "blue")
			//brdf = calc_phong_BRDF(albedo, albedo_spec, cos_theta);

		if(intersect_light(&shadow_ray_ip, shadow_ray))
		{
			float cos_light_angle = dot(shadow_ray_ip.normal, -shadow_ray.direction);
			float light_angle = shadow_ray_ip.area / n_samples * clamp(cos_light_angle, 0.0f, 1.0f) / pow(length(differance), 2) / (M_PI * 2);
			//float light_angle = shadow_ray_ip.area / n_samples * clamp(cos_light_angle, 0.0f, 1.0f) / pow(length(differance), 2) ;

			direct +=brdf *shadow_ray_ip.radiosity *cos_theta *light_angle;
		}
	}
	
	return direct;
}

Color Scene::global_diffuse(Ray r,int render_mode,InterPt ip,int iteration)
{
	Color global;
	static const int n_samples = 1;
	int temp_cotr = 0;
	for (int i = 0; i < n_samples; ++i)
	{
		// helper is just a random vector and can not possibly be
		// a zero vector since ip.normal is normalized
		vec3 helper = ip.normal + vec3(1,1,1);
		vec3 tangent = normalize(cross(ip.normal, helper));

		// rand1 is a random number from the cosine estimator
		float rand1 = (*dis_)(*gen_);
		float rand2 = (*dis_)(*gen_);

		// Uniform distribution over a hemisphere
		float inclination = acos(sqrt(rand1));//acos(1 - rand1);//acos(1 -  2 * (*dis_)(*gen_));
		float azimuth = 2 * M_PI * rand2;
		// Change the actual vector
		vec3 random_direction = ip.normal;
		random_direction = normalize(rotate(random_direction,inclination,tangent));
		random_direction = normalize(rotate(random_direction,azimuth,ip.normal));

		float cos_angle = dot(random_direction, ip.normal);
		float g = cos_angle / M_PI;

		Color brdf;
		brdf = calc_lambertian_BRDF(ip.material.color_diffuse * ip.material.reflectance * (1 - ip.material.specular_reflectance));

		r.direction = random_direction;
		r.radiance *= M_PI * brdf; // Importance sampling
		//r.radiance *= brdf; // Importance sampling
		//global += trace_ray(r, render_mode, iteration + 1) * M_PI * brdf;
		global += trace_ray(r, render_mode, iteration + 1) * r.radiance;
	}
	//cout<<"temp_cotr = "<<temp_cotr<<endl;
	return global / n_samples;
}

Color Scene::trace_mirror_reflection(Ray r,int render_mode,InterPt ip,int iteration)
{
	r.has_intersected = true;
	Color specular = Color();
	r.direction = reflect(r.direction, ip.normal);
	Color brdf = calc_ideal_BRDF(ip.material.color_specular * ip.material.reflectance * ip.material.specular_reflectance);
	r.radiance *= brdf;
	// Recursively trace the reflected ray
	specular += trace_ray(r, render_mode, iteration + 1) * brdf;
	return specular;
}

Color Scene::trace_refracted_ray(Ray r,int render_mode,InterPt ip,int iteration,vec3 offset,bool inside)
{
	Ray recursive_ray = r;
	recursive_ray.has_intersected = true;

	vec3 normal = inside ? -ip.normal : ip.normal;
	vec3 ideal_refraction = refract(r.direction,normal,r.material.refraction_index / ip.material.refraction_index);
	vec3 ideal_reflection = reflect(r.direction, ip.normal);
	if (ideal_refraction != vec3(0))
	{
		// https://en.wikipedia.org/wiki/Schlick%27s_approximation
		float n1 = r.material.refraction_index;
		float n2 = ip.material.refraction_index;
		float R_0 = pow((n1 - n2)/(n1 + n2), 2);
		float R = R_0 + (1 - R_0) * pow(1 - dot(normal, -r.direction),5);

		Ray recursive_ray_reflected = recursive_ray;
		Ray recursive_ray_refracted = recursive_ray;

		if (inside)
			offset = -offset;
		
		// Reflected ray
		// Change the material the ray is travelling in
		// https://computergraphics.stackexchange.com/questions/2482/choosing-reflection-or-refraction-in-path-tracing
		// Fresnel function https://en.wikipedia.org/wiki/Fresnel_equations
		recursive_ray_reflected.material = Material::air();
		recursive_ray_reflected.origin = r.origin + ip.t * r.direction +offset;
		// Refracted ray
		// Change the material the ray is travelling in
		recursive_ray_refracted.material = ip.material;
		recursive_ray_refracted.origin = r.origin + ip.t * r.direction -offset;
		
		Color to_return;
		recursive_ray_reflected.direction = ideal_reflection;
		recursive_ray_refracted.direction = ideal_refraction;

		Color brdf_specular(ip.material.color_specular * ip.material.reflectance * ip.material.specular_reflectance * R);
		Color brdf_refractive(ip.material.color_diffuse * ip.material.reflectance * ip.material.specular_reflectance * (1 - R));

		recursive_ray_reflected.radiance *= brdf_specular;
		recursive_ray_refracted.radiance *= brdf_refractive;

		// Recursively trace the refracted rays
		Color reflection = trace_ray(recursive_ray_reflected, render_mode, iteration + 1) * brdf_specular;
		Color refraction = trace_ray(recursive_ray_refracted, render_mode, iteration + 1) * brdf_refractive;
		return reflection + refraction;
	}
	else
	{ 
		if (inside)
			recursive_ray.origin = r.origin + ip.t * r.direction - offset;
		else
			recursive_ray.origin = r.origin + ip.t * r.direction + offset;

		Color brdf_specular(ip.material.color_specular * ip.material.reflectance * ip.material.specular_reflectance);
		recursive_ray.direction = ideal_reflection;
		recursive_ray.radiance *= brdf_specular;
		// Recursively trace the reflected ray
		return trace_ray(recursive_ray, render_mode, iteration + 1) * brdf_specular;
	}
}

Color Scene::trace_ray(Ray r, int render_mode, int iteration)
{
	InterPt ip;
	LightInterPt light_ip;

	if (intersect_light(&light_ip, r)) // Ray hit light source
		if (render_mode == SPECULAR)
			return light_ip.radiosity / (M_PI * 2);
			//return light_ip.radiosity;
		else
			return Color();

	else if (intersect(&ip, r))
	{ 
		float random = (*dis_)(*gen_);
		float esc_cond = iteration == 0 ? 1.0 : 0.8;
		if (random > esc_cond || iteration > 25)
			return Color(); // stop tracing
		vec3 offset = ip.normal * 0.00001f;
		bool inside = false;
		if (dot(ip.normal, r.direction) > 0) // The ray is inside an object
			inside = true;

		Color sum, reflection, diffuse;
		if (!ip.material.transmissivity) // glass 1 mirror 0
		{ //mirror or solid diffuse
			Ray recursive_ray = r;
			recursive_ray.origin = r.origin + ip.t * r.direction + (inside ? -offset : offset);
			reflection = Color(); 
			if(ip.material.specular_reflectance)
				reflection = trace_mirror_reflection(recursive_ray,render_mode,ip,iteration); //mirror - calculate specular reflection
			diffuse = Color();
			if (render_mode == MONTE_CARLO && ip.material.specular_reflectance == 0)
				diffuse = calc_diffuse(recursive_ray,render_mode,ip,iteration); // solid - calculate diffuse color
			//sum += (specular + diffuse) *(1 - ip.material.transmissivity); // if it's glass then -> transparent
			sum += (reflection + diffuse);
		} else { // glass
			Color transmitted = trace_refracted_ray(r, render_mode, ip, iteration, offset, inside);
			sum += transmitted * ip.material.transmissivity;
		}
		return sum / esc_cond;
	}
	return Color();
}
