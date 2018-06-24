#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <time.h>
#include <omp.h>

using namespace std;
#define Float double

struct Random {
	std::mt19937 engine;
	std::uniform_real_distribution<Float> dist;
	Random() {};
	Random(int seed) { engine.seed(seed); dist.reset(); }
	Float next() { return dist(engine); }
};

struct V
{
	Float x;
	Float y;
	Float z;
	V(Float v = 0)
		:V(v, v, v) {}
	V(Float x, Float y, Float z)
		:x(x), y(y), z(z) {}
	Float operator[](int i) const {
		return (&x)[i];
	}
};

int tonemap(Float v)
{// v^1/2.2 clamped to 0,255
	return std::min(std::max(int(std::pow(v, 1 / 2.2) * 255), 0), 255);
}

V operator+(V a, V b) {
	return V(a.x + b.x, a.y + b.y, a.z + b.z);
}
V operator-(V a, V b) {
	return V(a.x - b.x, a.y - b.y, a.z - b.z);
}
V operator*(V a, V b) {
	return V(a.x * b.x, a.y * b.y, a.z * b.z);
}
V operator/(V a, V b) {
	return V(a.x / b.x, a.y / b.y, a.z / b.z);
}
V operator-(V v) {
	return V(-v.x, -v.y, -v.z);
}

double dot(V a, V b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
V cross(V a, V b) {
	return V(a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}
V normalize(V v) {
	return v / sqrt(dot(v, v));
}

struct Ray
{
	V o;// origin
	V d;// direction (must be normalized)
};

struct Hit
{
	Float t;//distance from origin of the ray
	V p;//position
	V n;//normal (must be normalized)
	unsigned long long geomID;
};

struct Sphere
{
	V p;		//position:
	Float r;	//radius
	V R;		//reflectance
	V Le;		//illuminance
	Hit intersect(const Ray& ray, Float tmin, Float tmax) const
	{
		const V op = p - ray.o;
		const Float b = dot(op, ray.d);
		const Float det = b * b - dot(op, op) + r * r;
		if (det < 0) { return Hit{ -1.f,{},{},0 }; }
		const Float t1 = b - sqrt(det);
		if (tmin < t1&&t1 < tmax) {
			return Hit{ t1, {}, {}, 0 };
		}

		const Float t2 = b + sqrt(det);
		if (tmin < t2&&t2 < tmax) {
			return Hit{ t2, {}, {}, 0 };
		}
		return Hit{ -1.f,{},{},0 };
	}
};

struct Scene
{
	std::vector<Sphere> spheres{
		/*  { V(0,-1,4)       , 1., V(0.7,0,0)  },
		  { V(0,1,4)       , 0.71, V(0.7,0.5,0)  },*/
		  { V(1e5 + 1,40.8,81.6)  , 1e5 , V(.75,.25,.25) },
		  { V(-1e5 + 99,40.8,81.6), 1e5 , V(.25,.25,.75) },
		  { V(50,40.8,1e5)      , 1e5 , V(.75) },
		  { V(50,1e5,81.6)      , 1e5 , V(.75) },
		  { V(50,-1e5 + 81.6,81.6), 1e5 , V(.75) },
		  { V(27,16.5,47)       , 16.5, V(.999)  },
		  { V(73,16.5,78)       , 16.5, V(.999)  },
		  { V(50,681.6 - .27,81.6), 600 , V(), V(12) },
	};
	Hit intersect(const Ray& ray, Float tmin, Float tmax) const
	{
		Hit minh;
		for (unsigned long long i = 0; i < spheres.size(); i++)
		{
			const Sphere& sphere = spheres.at(i);
			const Hit h = sphere.intersect(ray, tmin, tmax);
			if (h.t > 0.) {
				minh = h;
				minh.geomID = i;
				tmax = minh.t;
			}
		}

		if (minh.t > 0.) {
			const Sphere& s = spheres.at(minh.geomID);
			minh.p = ray.o + ray.d * minh.t;
			minh.n = (minh.p - s.p) / s.r;//should be normalized thanks to s.r but floating point errors can occur
		}
		else
			minh.geomID = -1;
		return minh;
	}
};

int main()
{
	//Camera
	const V c0 = V(50, 52, 295.6);//origin of rays
	const V poi = c0 + V(0, -0.042612, -1);//point of interest
	const V up(0, 1, 0);//vector up for the camera

	const V dir = normalize(c0 - poi);
	const V u = normalize(cross(up, dir));//for x vector
	const V v = cross(dir, u);//dir and u are normalized so no need to normalize here

	int w = 1200, h = 800, spp=8;
	Float fovx = 30.*M_PI / 180.;
	Float dd = (Float)(tan(fovx*0.5));
	Float ddx = (Float)(tan(fovx*0.5))*2. / w;//size of a pixel
	Float ar = (h > w) ? (Float)w / h : (Float)h / w;
	//other params


	Scene s;

	std::vector<V> I(w*h);

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < w*h; i++)
	{
		thread_local Random rng(42 + omp_get_thread_num());
		for (int j = 0; j < spp; j++) {
			const int x = i % w;
			const int y = i / w;

			Ray ray;
			//Perspective
			//ray.o = V(0,0,0);
			//ray.d = normalize(V(x*dd - dd*(Float)w/2. + dd/2., -y*dd + dd*(Float)h/2. - dd/2, 1.));//center of the pixel

			ray.o = c0;
			V vw = ((x+rng.next())*ddx - ddx * w*0.5 + ddx * 0.5)*u + (ddx*h*0.5 - ddx * 0.5 - ddx * (y+rng.next()))*v - ar * dir;
			ray.d = normalize(vw);
			//V w( (2.*(Float)x / w - 1.)*dd*ar, dd*(1.-2.*(Float)y / h), -1.);
			//w = normalize(w);
			//ray.d = w.x*u + w.y*v + w.z*dir;

			//Orthographique
			// ray.o = V( 2.*((Float)x/w)-1,2.*((Float)y/h)-1,5);
			// ray.d = V(0,0,-1);

			Hit hit = s.intersect(ray, 0., 1e+10);
			if (hit.geomID == -1)
			{
				I[y*w + x] = I[y*w + x]+ V(0)/spp;
			}
			else
			{
				V c;
				//c=hit.n;
				Float cosTerm = dot(hit.n, -ray.d);
				c = s.spheres.at(hit.geomID).R *cosTerm;
				I[y*w + x] = I[y*w +x] + c/spp;
			}
		}

	}

	ofstream ofs;
	ofs.open("image.ppm", ios_base::trunc | ios_base::binary);//delete the image at the opening
	ofs << "P3 " << w << " " << h << " 255\n";

	for (const auto& i : I)
	{
		ofs << tonemap(fabs(i.x)) << " " << tonemap(fabs(i.y)) << " " << tonemap(fabs(i.z)) << "\n";
	}
	ofs.close();

	return 0;
}

