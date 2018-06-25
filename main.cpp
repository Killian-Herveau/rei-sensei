#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <time.h>
#include <omp.h>
#include <tuple>

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
V sample_sphere(Float u, Float v)
{
	V r;
	r.z = 1. - 2.* v;
	Float sinT = sqrt(max(0., 1. - r.z*r.z));
	Float phi = 2.*Float(M_PI)*u;
	r.x = sinT * cos(phi);
	r.y = sinT * sin(phi);
	return r;
}

std::tuple<V, V> tangentSpace(const V& n) {
	const double s = std::copysign(1, n.z);
	const double a = -1 / (s + n.z);
	const double b = n.x*n.y*a;
	return {
		V(1 + s * n.x*n.x*a,s*b,-s * n.x),
		V(b,s + n.y*n.y*a,-n.y)
	};
}

V sample_hemisphere(Float u, Float v, V n)
{
	V r = sample_sphere(u, v);
	return normalize((dot(r, n) > 0.) ? r : -r);
}

V sample_hemisphere2(Float r1, Float r2, V n)
{
	const tuple<V, V>& tuple = tangentSpace(n);
	const double r = sqrt(r1);
	const double t = 2 * M_PI*r2;
	const double x = r * cos(t);
	const double y = r * sin(t);
	const auto d = V(x, y, std::sqrt(std::max(.0, 1 - x * x - y * y)));
	// Convert to world coordinates
	return tuple._Myfirst._Val * d.x + tuple._Get_rest()._Myfirst._Val * d.y + n * d.z;

}

struct Ray
{
	V o;// origin
	V d;// direction (must be normalized)
};
enum class SurfaceType {
	Diffuse,
	Mirror,
	Fresnel,
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
	SurfaceType s;
	V R;		//reflectance
	V Le;		//illuminance
	Float ior = 1.5168;
	Hit intersect(const Ray& ray, Float tmin, Float tmax) const
	{
		const V op = p - ray.o;
		const Float b = dot(op, ray.d);
		const Float det = b * b - dot(op, op) + r * r;
		if (det < 0) { return Hit{ -1.,{},{},0 }; }
		const Float t1 = b - sqrt(det);
		if (tmin < t1&&t1 < tmax) {
			return Hit{ t1, {}, {}, 0 };
		}

		const Float t2 = b + sqrt(det);
		if (tmin < t2&&t2 < tmax) {
			return Hit{ t2, {}, {}, 0 };
		}
		return Hit{ -1.,{},{},0 };
	}
};

struct Scene
{
	std::vector<Sphere> spheres{
		/*  { V(0,-1,4)       , 1., V(0.7,0,0)  },
		  { V(0,1,4)       , 0.71, V(0.7,0.5,0)  },*/
		  { V(1e5 + 1,40.8,81.6)	, 1e5 ,SurfaceType::Diffuse, V(.75,.25,.25) },
		  { V(-1e5 + 99,40.8,81.6)	, 1e5 ,SurfaceType::Diffuse, V(.25,.25,.75) },
		  { V(50,40.8,1e5)			, 1e5 ,SurfaceType::Diffuse, V(.75) },
		  { V(50,1e5,81.6)			, 1e5 ,SurfaceType::Diffuse, V(.75) },
		  { V(50,-1e5 + 81.6,81.6)	, 1e5 ,SurfaceType::Diffuse, V(.75) },
		  { V(27,16.5,47)			, 16.5,SurfaceType::Mirror,  V(.999)  },
		  
		  { V(73,16.5,78)			, 16.5,SurfaceType::Fresnel, V(.999)},
		  { V(50,93,50)			    , 20.5 ,SurfaceType::Diffuse, V(.999), V(12) },
		  //{ V(50,681.6 - .27,81.6)	, 600 ,SurfaceType::Diffuse, V(), V(12) },
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
	const V c0 = V(50, 52, 295.6);			//origin of rays
	const V poi = c0 + V(0, -0.042612, -1);	//point of interest
	const V up(0, 1, 0);					//vector up for the camera

	const V dir = normalize(c0 - poi);
	const V u = normalize(cross(up, dir));	//for x vector
	const V v = cross(dir, u);				//dir and u are normalized so no need to normalize here
	clock_t cbegin = clock();
	int w = 1200, h = 800, spp = 2000;
	Float fovx = 30.*M_PI / 180.;
	Float dd = (Float)(tan(fovx*0.5));
	Float ddx = (Float)(tan(fovx*0.5))*2. / w;//size of a pixel
	Float ar = (h > w) ? (Float)w / h : (Float)h / w;
	//other params

	Scene s;
	std::vector<V> I; I.resize(w*h, V(-0.1415));

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < w*h; i++)
	{
		I[i] = 0;//because it is initialized at -0.1415 for counting %
		thread_local Random rng(42 + omp_get_thread_num());
		for (int j = 0; j < spp; j++) {
			const int x = i % w;
			const int y = i / w;

			Ray ray;
			//Perspective
			//ray.o = V(0,0,0);
			//ray.d = normalize(V(x*dd - dd*(Float)w/2. + dd/2., -y*dd + dd*(Float)h/2. - dd/2, 1.));
			//center of the pixel

			ray.o = c0;
			V vw = ((x + rng.next())*ddx - ddx * w*0.5 + ddx * 0.5)*u +  // u direction
				(ddx*h*0.5 - ddx * 0.5 - ddx * (y + rng.next()))*v + // v direction
				-ar * dir;											 // dir direction
			ray.d = normalize(vw);

			//V w( (2.*(Float)x / w - 1.)*dd*ar, dd*(1.-2.*(Float)y / h), -1.);
			//w = normalize(w);
			//ray.d = w.x*u + w.y*v + w.z*dir;
			//Orthographic
			// ray.o = V( 2.*((Float)x/w)-1,2.*((Float)y/h)-1,5);
			// ray.d = V(0,0,-1);

			V L(0), th(1);//Luminance, throughput
			for (int depth = 0; depth < 10; depth++)
			{
				Hit hit = s.intersect(ray, 1e-4, 1e+10);
				if (hit.geomID == -1)
				{//background color

					I[y*w + x] = I[y*w + x] + V(0) / spp;
					break;//will not be reflected or anything. V(0) is background color
				}
				
				Sphere& sph = s.spheres.at(hit.geomID);
				V c;
				L = L + th * sph.Le;// *cosTerm*M_1_PI;

				ray.o = hit.p;
				if (sph.s == SurfaceType::Diffuse)
				{
					if (dot(hit.n, -ray.d) < 0)
						hit.n = -hit.n;
					ray.d = sample_hemisphere(rng.next(), rng.next(), hit.n);
          th = th * sph.R * dot(hit.n, ray.d) * 2; //updating throughput
				}
				else if (sph.s == SurfaceType::Mirror)//reflect
				{
					if (dot(hit.n, -ray.d) < 0)
						hit.n = -hit.n;
					ray.d = 2.*dot(-ray.d, hit.n) * hit.n + ray.d;
          th = th * sph.R; //updating throughput
				}
				else if (sph.s == SurfaceType::Fresnel)
				{
					const V wi = -ray.d;
					const bool into = dot(wi, hit.n) > 0;
					const V n = (into) ? hit.n : -hit.n;
					const Float eta = (into) ? 1. / sph.ior : sph.ior;

					const Float t1 = dot(wi, n);
					const Float t2 = 1. - eta * eta*(1. - t1 * t1);

					if (t2 < 0)//total internal reflection
						ray.d = 2.*dot(wi, hit.n)*hit.n - wi;
					else
					{
						const V wt = eta * (n*t1 - wi) - n * sqrt(t2);
						const Float c = (into) ? dot(wi, hit.n) : dot(wt, hit.n);
						const Float r = (1. - sph.ior) / (1. + sph.ior);
						const Float Fr = r * r + (1. - r * r)*pow(1 - c, 5); // Schlick
						// Fr term tells us if it's more likely to have refraction or reflection
						// and we choose randomly one path or the other
						ray.d = (rng.next() < Fr) ? 2.*dot(wi, hit.n)*hit.n - wi : wt;
					}
          th = th * sph.R; //updating throughput
				}

				if (max({ th.x,th.y,th.z }) == 0.)
					break; // if throughput is too low
			}
			I[i] = I[i] + L / spp;

		}

		if (rng.next() > 0.99997)
		{
			int count = 0;
			for (const V& v : I)
				if (v.x == -0.1415)
					count++;
			cout << "completion : " << (w*h - count)*100. / (w*h) << "%\n";
		}
	}
	std::cout << "time : " << double(clock() - cbegin) / CLOCKS_PER_SEC;
	ofstream ofs;
	ofs.open("image.ppm", ios_base::trunc);//delete the image at the opening
	ofs << "P3 " << w << " " << h << " 255\n";

	for (const auto& i : I)
	{
		ofs << tonemap(abs(i.x)) << " " << tonemap(abs(i.y)) << " " << tonemap(abs(i.z)) << "\n";
	}
	ofs.close();

	return 0;
}

