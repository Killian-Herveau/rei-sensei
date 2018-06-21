#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <omp.h>
using namespace std;
#define Float float

struct V
{
    Float x;
    Float y;
    Float z;
    V(Float v=0)
        :V(v,v,v){}
    V(Float x,Float y,Float z)
        :x(x),y(y),z(z) {}
    Float operator[](int i) const {
        return (&x)[i];
    }
};

int tonemap(Float v)
{// v^1/2.2 clamped to 0,255
    return std::min(std::max(int(std::pow(v,1/2.2)*255),0),255);
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
    V o;
    V d;
};

struct Hit
{
    Float t;//distance from origin of the ray
    V p;//position
    V n;//normal
    unsigned long long geomID;
};

struct Sphere
{
    V p;//position:
    Float r;//radius
    V R;
    V Le;
    Hit intersect(const Ray& ray, Float tmin, Float tmax) const
    {
        const V op = p - ray.o;
        const Float b = dot(op, ray.d);
        const Float det = b * b - dot(op, op) + r*r;
        if (det < 0) { return {}; }
        const Float t1 = b - sqrt(det);
        if (tmin<t1&&t1<tmax) {
            return Hit{t1, {}, {}, 0}; }

        const Float t2 = b + sqrt(det);
        if (tmin<t2&&t2<tmax) {
            return Hit{t2, {}, {}, 0}; }
        return Hit{-1.f,{},{},0};
    }
};

struct Scene
{
    std::vector<Sphere> spheres{
        { V(0,0,3)       , 1, V(0.7,0,0)  },
        /* { V(1e5+1,40.8,81.6)  , 1e5 , V(.75,.25,.25) },
        { V(-1e5+99,40.8,81.6), 1e5 , V(.25,.25,.75) },
        { V(50,40.8,1e5)      , 1e5 , V(.75) },
        { V(50,1e5,81.6)      , 1e5 , V(.75) },
        { V(50,-1e5+81.6,81.6), 1e5 , V(.75) },
        { V(27,16.5,47)       , 16.5, V(.999)  },
        { V(73,16.5,78)       , 16.5, V(.999)  },
        { V(50,681.6-.27,81.6), 600 , V(), V(12) },*/
    };
    Hit intersect(const Ray& ray, Float tmin, Float tmax) const
    {
        Hit minh;
        for(unsigned long long i = 0 ; i < spheres.size() ; i++)
        {
            const Sphere& sphere=spheres.at(i);
            const Hit h = sphere.intersect(ray, tmin, tmax);
            if (h.t<0.) { continue; }
            minh = h;
            tmax = minh.t;
            minh.geomID=i;
        }

        if (minh.t>0.) {
            const Sphere& s = spheres.at(minh.geomID);
            minh.p = ray.o + ray.d * minh.t;
            minh.n = (minh.p - s.p) / s.r;
        }
        else
            minh.geomID = -1;
        return minh;
    }
};

int main()
{
    //Camera
    const V c0 = V(0,0,0);//origin of rays
    const V poi = V(0,0,1);//point of interest
    const V up(0,1,0);//vector up for the camera

    const V dir = normalize(poi-c0);
    const V u = normalize(cross(up,dir));//for x vector
    const V v = cross(dir,u);//dir and u are normalized so no need to normalize here

    int w=1280, h=720;
    float fovx = 80.*M_PI/180.;
    Float dd = 2.*(Float)(tan(fovx*0.5))/w;

    //other params
    ofstream ofs;
    ofs.open("image.ppm",ios_base::trunc | ios_base::binary);//delete the image at the opening



    ofs<<"P3 "<<w<<" "<<h<<" 255\n";

    Scene s;
    for(int i = 0; i< w*h ; i++)
    {
        const int x = i % w;
        const int y = i / w;

        Ray ray;
        //Perspective
        //ray.o = V(0,0,0);
        //ray.d = normalize(V(x*dd - dd*(Float)w/2. + dd/2., -y*dd + dd*(Float)h/2. - dd/2, 1.));//center of the pixel

        ray.o = c0;
        V test = c0+dir+(x*dd-dd*w*0.5+dd*0.5)*u + (dd*h*0.5 - dd*0.5 - dd*y)*v;
        ray.d = normalize(test);

        //Orthographique
        // ray.o = V( 2.*((Float)x/w)-1,2.*((Float)y/h)-1,5);
        // ray.d = V(0,0,-1);

        Hit hit = s.intersect(ray,0.,1e+10);
        if(hit.geomID==-1)
            ofs<<"0 0 0\n";
        else
        {
            V c;
            c=hit.n;
            ofs<<tonemap(c.x) << " " << tonemap(c.y) << " " << tonemap(c.z) << "\n";
         }
        //ofs<<"150 0 0\n";
    }
    ofs.close();

    return 0;
}

