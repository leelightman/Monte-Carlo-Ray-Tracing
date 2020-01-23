#include <string>
#include <iostream>
#include <Eigen/Dense>

#include "Color.h"
using namespace glm;
using namespace std;
using namespace Eigen;

Color::Color()
{
	for (int i = 0; i < NUM_CHANNELS; ++i) //NUM_CHANNELS is 3
	{
		data[i] = 0;
	}
}

float Color::norm() const //average
{
	float sum = 0;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		sum += data[i];
	}
	return sum / NUM_CHANNELS;
}

ostream& operator<<(ostream& os, const Color& sd)
{
	os << "[ ";
	for (int i = 0; i < sd.NUM_CHANNELS - 1; ++i)
	{
		os << sd.data[i] << ", ";
	}
	os << sd.data[sd.NUM_CHANNELS - 1] << "]";
	return os;
}

Color operator*(float f, const Color& sd)
{
	return sd * f;
}

float& Color::operator[](const int i)
{
	return data[i];
}

Color Color::operator+(const Color& sd) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = data[i] + sd.data[i];
	}
	return to_return;
}

Color Color::operator-(const Color& sd) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = data[i] - sd.data[i];
	}
	return to_return;
}

Color Color::operator^(const float& f) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = pow(data[i], f);
	}
	return to_return;
}

Color Color::operator/(const float& f) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = data[i] / f;
	}
	return to_return;
}

Color Color::operator*(const float& f) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = data[i] * f;
	}
	return to_return;
}

Color Color::operator*(const Color& sd) const
{
	Color to_return;
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		to_return.data[i] = data[i] * sd.data[i];
	}
	return to_return;
}


Color Color::operator+=(const Color& sd)
{
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		data[i] = data[i] + sd.data[i];
	}
	return *this;
}

Color Color::operator-=(const Color& sd)
{
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		data[i] = data[i] - sd.data[i];
	}
	return *this;
}

Color Color::operator*=(const Color& sd)
{
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		data[i] = data[i] * sd.data[i];
	}
	return *this;
}


Color Color::operator/=(const float& f)
{
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		data[i] = data[i] / f;
	}
	return *this;
}

Color Color::operator*=(const float& f)
{
	for (int i = 0; i < NUM_CHANNELS; ++i)
	{
		data[i] = data[i] * f;
	}
	return *this;
}

Color Color::white()
{
	Color white;
	white[0] = 1;
	white[1] = 1;
	white[2] = 1;
	return white;
}

Material Material::air()
{
	Material air;
	Color color_diffuse;
	Color color_specular;
	air.color_diffuse = color_diffuse;
	air.color_specular = color_specular;
	air.reflectance = 0;
	air.specular_reflectance = 0;
	air.transmissivity = 1;
	air.refraction_index = 1;
	return air;
}

//https://www.cs.cmu.edu/afs/cs/academic/class/15462-f09/www/lec/lec8.pdf

Color calc_lambertian_BRDF(Color albedo)
{
	return albedo / M_PI;
}
Color calc_ideal_BRDF(Color albedo)
{
	return albedo;
}

Color calc_phong_BRDF(Color albedo, Color albedo_spec, float cos_theta)
{
	const float PHONG_COEFFICIENT = 1000;
	return albedo/M_PI+albedo_spec*(PHONG_COEFFICIENT+2)/M_PI*pow(cos_theta,PHONG_COEFFICIENT);
}