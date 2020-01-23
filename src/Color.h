#ifndef COLOR_H
#define COLOR_H

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <Eigen/Dense>

using namespace std;
using namespace glm;
using namespace Eigen;

class Color
{
public:
	Color();

	float norm() const;

	// Various operators for Color
	friend ostream& operator<<(
		ostream& os,
		const Color& sd);
	friend Color operator*(
		float f,
		const Color& sd);
	float& operator[](const int i);
	Color operator+(const Color& sd) const;
	Color operator-(const Color& sd) const;
	Color operator^(const float& f) const;
	Color operator/(const float& f) const;
	Color operator*(const float& f) const;
	Color operator*(const Color& sd) const;
	Color operator+=(const Color& sd);
	Color operator-=(const Color& sd);
	Color operator*=(const Color& sd);
	Color operator/=(const float& f);
	Color operator*=(const float& f);

	// Currently not used as wavelengths. We only care about three channels,
	// (hence three wavelengths) r, g, b. These would not correspond to real
	// wavelengths at the moments since r is the lowest and b is the highest.
	static const int NUM_CHANNELS = 3;
	// The data contains radiance values in all wave lengths.
	float data[NUM_CHANNELS];
	static Color white();
};

struct Material
{
	// The colors are reflectance distribution functions
	// All channels are in the interval [0 , 1]
	Color color_diffuse;
	Color color_specular;
	float reflectance; 
	float specular_reflectance; 
	float transmissivity;
	float refraction_index;
	string name;

	static Material air();
};

struct Ray
{
	vec3 origin;
	vec3 direction;
	Material material; // The material the ray is travelling in
	// The radiance variable is used when tracing rays from light source to photons.
	// When tracing from the camera, the radiance variable will be used for importance. 
	Color radiance; // [Watts / m^2 / steradian]
	bool has_intersected;  // This is used only when forward tracing ray
};

struct InterPt
{
	Material material; // Material of the object hit by the ray
	vec3 normal; // Normal of the surface hit by the ray
	float t; // The distance the ray travelled before intersecting
};

struct LightInterPt
{
	Color radiosity; // The radiosity of the light source [Watts/m^2]
	float area; // The area of the light source [m^2]
	vec3 normal; // Normal of the surface hit by the ray
	float t; // The distance the ray travelled before intersecting
};

Color calc_ideal_BRDF(Color albedo);
//Color calcLambertianBRDF(vec3 d1,vec3 d2,vec3 normal,Color albedo);
Color calc_lambertian_BRDF(Color albedo);
Color calc_phong_BRDF(Color albedo, Color albedo_spec, float cos_theta);

#endif // Camera_H