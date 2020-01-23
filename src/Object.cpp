#include "Object.h"
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
using namespace std;
using namespace glm;
using namespace Eigen;


struct Vertex{
    int index;
    float x;
    float y;
    float z;
};

struct Face{
    int v_0;
    int v_1;
    int v_2;
};

struct Model{
    Array<Vertex,502,1> V;
    Array<Face,1000,1> F;
};

struct Cube{
    Array<Vertex,8,1> V;
    Array<Face,10,1> F;
};

Model read_off_file(string filename){
    string readLine;
    ifstream inFile;
    Model model_instance;
    inFile.open(filename);
    if(!inFile){
        cerr << "Unable to open"<<filename<<endl;
        exit(1);   // call system to stop
    } else{
        getline(inFile, readLine);
        //if (readLine != "OFF"){
            //cout << "The file"<<filename<<" to read is not in OFF format." << endl;
        cout<<filename<<endl;
        cout<<readLine<<endl;
        //}
        getline(inFile, readLine);
        
        // Find out number of vertices
        int space1 = readLine.find(" ", 0);
        int vertex_number = atoi(readLine.substr(0,space1+1).c_str());
        
        // Find out number of faces
        int space2 = readLine.find(" ", space1+1);
        int face_number = atoi(readLine.substr(space1,space2+1).c_str());
        int space3 = 0;
        int space4 = 0;

        // Read vertices into array V
        Array<Vertex, 502,1> V;
        for (int n=0; n<vertex_number; n++){
            getline(inFile, readLine);
            V[n].index = n;
            space1 = readLine.find(" ", 0); 
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find("\n", space2+1);

            V[n].x = atof(readLine.substr(0,space1+1).c_str());
            V[n].y = atof(readLine.substr(space1,space2+1).c_str());
            V[n].z = atof(readLine.substr(space2).c_str());
        }

        // read faces into array F
        Array<Face, 1000,1> F;
        for (int n=0; n<face_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0);
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find(" ", space2+1);
            space4 = readLine.find(" ", space3+1);

            F[n].v_0 = atoi(readLine.substr(space1,space2+1 ).c_str());            
            F[n].v_1 = atoi(readLine.substr(space2,space3+1 ).c_str());            
            F[n].v_2 = atoi(readLine.substr(space3).c_str());   
        }
        model_instance.V = V;
        model_instance.F = F;
    }
    inFile.close();
    //cout<<"step"<<endl;
    return model_instance;
}

Object::Object(Material* material) : 
	material_(material)
{}

Material Object::material() const
{
	return material_ ? *material_ : Material();
}

//Sphere

Sphere::Sphere(vec3 position, float radius, Material* material) : 
	Object(material), POSITION_(position), RADIUS_(radius)
{}

bool Sphere::intersect(InterPt* id, Ray r) const
{
	// if to_square is negative we have imaginary solutions,
	// hence no intersection
	// p_half comes from the p-q formula (p/2)
	float p_half = dot((r.origin - POSITION_), r.direction);
	float to_square = pow(p_half, 2) + pow(RADIUS_, 2) - pow(length(r.origin - POSITION_), 2);
	float t; // parameter that tells us where on the ray the intersection is
	vec3 n; // normal of the intersection point on the surface
	if (to_square < 0)
	// No intersection points
		return false;
	else // if (to_square > 0) or (to_square == 0)
	// One or two intersection points, if two intersection points,
	// we choose the closest one that gives a positive t
	{
		t = -p_half - sqrt(to_square); // First the one on the front face
		if (t < 0) // if we are inside the sphere
		{
			// the intersection is on the inside of the sphere
			t = -p_half + sqrt(to_square);
		}
		n = r.origin + t*r.direction - POSITION_;
	}
	if (t >= 0) // t needs to be positive to travel forward on the ray
	{
		id->t = t;
		id->normal = normalize(n);
		id->material = material();
		return true;
	}
	return false;
}

vec3 Sphere::get_pt_on_surface(float u, float v) const
{
	// Uniform over a sphere
	float inclination = acos(1 - 2 * u);
	float azimuth = 2 * M_PI * v;

	vec3 random_direction = vec3(1,0,0);
	random_direction = normalize(rotate(random_direction,inclination,vec3(0,1,0)));
	random_direction = normalize(rotate(random_direction,azimuth,vec3(1,0,0)));

	return POSITION_ + random_direction * RADIUS_;
}

// Plane

Plane::Plane(vec3 p0, vec3 p1, vec3 p2, Material* material) : 
	Object(material),
	P0_(p0),
	P1_(p1),
	P2_(p2),
	NORMAL_(normalize(cross(p0 - p1, p0 - p2))),
	AREA_(length(cross(p0 - p1, p0 - p2)))
{}

bool Plane::intersect(InterPt* id, Ray r) const
{
	const double EPSILON = 0.00001;

	// Möller–Trumbore intersection algorithm
	vec3 e1, e2;  //Edge1, Edge2
	vec3 P, Q, T;
	float det, inv_det, u, v;
	float t;

	// Find vectors for two edges sharing P0_
	e1 = P1_ - P0_;
	e2 = P2_ - P0_;
	P = cross(r.direction, e2);
	det = dot(e1, P);
	// NOT CULLING
	if(det > -EPSILON && det < EPSILON) return false;
	inv_det = 1.0 / det;

	// calculate distance from P0_ to ray origin
	T = r.origin - P0_;
	Q = cross(T, e1);

	// Calculate u parameter and test bound
	u = dot(T, P) * inv_det;
	v = dot(r.direction, Q) * inv_det;

	// The intersection lies outside of the plane
	if(u < 0.f || u > 1.f || v < 0.f || v > 1.f) return false;

	t = dot(e2, Q) * inv_det;

	if(t > EPSILON ) {
		id->t = t;
		id->normal = normalize(cross(e1, e2));
		id->material = material();
		return true;
	}

	return false;
}

vec3 Plane::get_pt_on_surface(float u, float v) const
{
	vec3 v1 = P1_ - P0_;
	vec3 v2 = P2_ - P0_;
	return P0_ + u * v1 + v * v2;
}

float Plane::get_area() const
{
	return AREA_;
}

vec3 Plane::get_normal() const
{
	return NORMAL_;
}

vec3 Plane::get_first_tangent() const
{
	return normalize(P1_ - P0_);
}

// Mesh
// load cube

bool read_cube(string filename, vector<vec3>& positions, vector<unsigned int>& faces){
    string readLine;
    ifstream inFile;
    inFile.open(filename);
    if(!inFile){
        cerr << "Unable to open"<<filename<<endl;
        exit(1);   // call system to stop
    } else{
        getline(inFile, readLine);
        //if (readLine != "OFF"){
            //cout << "The file"<<filename<<" to read is not in OFF format." << endl;
        //cout<<filename<<endl;
        //cout<<readLine<<endl;
        //}
        getline(inFile, readLine);
        
        // Find out number of vertices
        int space1 = readLine.find(" ", 0);
        int vertex_number = atoi(readLine.substr(0,space1+1).c_str());
        
        // Find out number of faces
        int space2 = readLine.find(" ", space1+1);
        int face_number = atoi(readLine.substr(space1,space2+1).c_str());
        int space3 = 0;
        int space4 = 0;

        // Read vertices into array V
        vec3 vertex;
        for (int n=0; n<vertex_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0); 
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find("\n", space2+1);

            vertex.x = atof(readLine.substr(0,space1+1).c_str());
            vertex.y = atof(readLine.substr(space1,space2+1).c_str());
            vertex.z = atof(readLine.substr(space2).c_str());
			positions.push_back(vertex);
        }

        // read faces into array F
        //Array<Face, 10,1> F;
        for (int n=0; n<face_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0);
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find(" ", space2+1);
            space4 = readLine.find(" ", space3+1);

            faces.push_back(atoi(readLine.substr(space1,space2+1 ).c_str()));            
            faces.push_back(atoi(readLine.substr(space2,space3+1 ).c_str()));            
            faces.push_back(atoi(readLine.substr(space3).c_str()));   
        }
    }
    inFile.close();
    return true;
}

Mesh::Mesh(mat4 transform, string filename, Material * material) :
	Object(material)
{
	transform_ = transform;

	vector<vec3> tmp_positions;
	vector<unsigned int> indices_;

	if(!read_cube(filename, tmp_positions, indices_))
		exit(2);
	for (int i = 0; i < tmp_positions.size(); ++i)
	{
		tmp_positions[i] = vec3(transform_ * vec4(tmp_positions[i], 1));
	}
	this->positions_ = tmp_positions;
	this->indices_ = indices_;
}

bool Mesh::intersect(InterPt* id, Ray r) const
{	
	//double t_local = 1000;
	for (int i = 0; i < this->indices_.size(); i=i+3)
	{
		vec3 v0 = this->positions_[i];
		vec3 v1 = this->positions_[i+1];
		vec3 v2 = this->positions_[i+2];
		vec3 edge1, edge2, h, s, q;
		double t = 0;
		const double EPSILON = 0.00001f;
		double a,f,u,v;
		edge1 = v1 - v0;
		edge2 = v2 - v0;
		h = cross(r.direction, edge2);
		a = dot(edge1,h);
		if (a > -EPSILON && a < EPSILON)
			return false;  
		f = 1.0/a;
		s = r.origin - v0;
		u = f * dot(s,h);
		if (u < 0.0 || u > 1.0)
			return false;
		q = cross(s,edge1);
		v = f * dot(r.direction,q);
		//if (v < 0.0 || u + v > 1.0)
			//return false;
		t = f * dot(edge2,q);
		//if (t > EPSILON && t < t_local){
		if (t > EPSILON){
			id->normal = r.origin + r.direction * t;
			id->t = t;
			id->material = material();
			//t_local = t;
			return true;
		} else
        	return false;
	}
	/*
	cout<<"tlocal: "<<t_local<<endl;
	if (t_local == 1000){
		return false;
	} else{
		return true;
	}*/
}


// Light
Light::Light(vec3 p0,vec3 p1,vec3 p2,float flux,Color color) :
	emitter_(p0, p1, p2, NULL),
	radiosity(flux / emitter_.get_area() * color)
{}

bool Light::intersect(LightInterPt* light_id, Ray r)
{
	InterPt id;
	if(emitter_.intersect(&id, r)) //intersect with the plane
	{
		light_id->normal = id.normal;
		light_id->t = id.t;
		light_id->radiosity = radiosity;
		light_id->area = get_area();
		return true;
	}
	else
		return false;
}

vec3 Light::get_pt_on_surface(float u, float v)
{
	return emitter_.get_pt_on_surface(u, v);
}

float Light::get_area() const
{
	return emitter_.get_area();
}

vec3 Light::get_normal() const
{
	return emitter_.get_normal();
}

Ray Light::shoot_light_ray()
{
	// Move random code out later
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(0, 1);

	Ray r;
	r.origin = get_pt_on_surface(dis(gen), dis(gen));

	// Get a uniformly distributed vector
	vec3 normal = emitter_.get_normal();
	vec3 tangent = emitter_.get_first_tangent();
	// rand1 is a random number from the cosine estimator
	float rand1 = dis(gen);
	float rand2 = dis(gen);

	// Uniform distribution
	float inclination = acos(sqrt(rand1));//acos(1 - rand1);//acos(1 -  2 * (*dis_)(*gen_));
	float azimuth = 2 * M_PI * rand2;
	// Change the actual vector
	vec3 random_direction = normal;
	random_direction = normalize(rotate(random_direction,inclination,tangent));
	random_direction = normalize(rotate(random_direction,azimuth,normal));

	r.direction = random_direction;
	r.material = Material::air();
	
	return r;
}