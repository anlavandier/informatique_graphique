#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <random>
#include <omp.h>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <chrono>
#include <array>
#include <memory>


static inline double sqr(double x) { return x * x; }

std::default_random_engine generator[64];
std::uniform_real_distribution<double> uniform(0.0,1.0);


// ====================================================================================================================
// Class Vector
// ====================================================================================================================
class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

	void normalize() {
		double n = 1 / sqrt(this->norm2());
		coord[0] *= n;
		coord[1] *= n;
		coord[2] *= n;
	}

	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, double b) {
	double over_b = 1/b;
	return Vector(a[0] * over_b, a[1] * over_b, a[2] * over_b);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator/(const Vector& a, const Vector& b) {
	return Vector(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}
Vector abs(const Vector& a) {
	return Vector(std::abs(a[0]), std::abs(a[1]), std::abs(a[2]));
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

std::ostream& operator<< (std::ostream& os, Vector &a) {
	os << a[0] << ' ' << a[1] << ' ' << a[2];
	return os;
}

Vector random_cos(const Vector& normal) {
    int thread_id = omp_get_thread_num();
    double r1 = uniform(generator[thread_id]);
    double r2 = uniform(generator[thread_id]);
    double r = sqrt(1 - r2);
    double x = cos(2. * M_PI * r1) * r;
    double y = sin(2. * M_PI * r1) * r;
    double z = sqrt(r2);

    Vector tangent1;
    if ((abs(normal[0]) <= abs(normal[1])) && (abs(normal[0]) <= abs(normal[2]))) {
        tangent1 = Vector(0, -normal[2], normal[1]);
    } else if ((abs(normal[1]) <= abs(normal[2])) && (abs(normal[1]) <= abs(normal[0]))) {
        tangent1 = Vector(-normal[2], 0, normal[0]);
    } else {
        tangent1 = Vector(-normal[1], normal[0], 0);
    }

    tangent1.normalize();
    Vector tangent2 = cross(normal, tangent1);
    return (z * normal + x * tangent1 + y * tangent2);
}


// ====================================================================================================================
// Class Ray
// ====================================================================================================================
class Ray {
public:
	explicit Ray(const Vector& o, const Vector& u) {
		origin = o;
		direction = u / sqrt(u.norm2());
	}

	Vector point_along(double t) {
		return origin + t * direction;
	}

	Vector origin;
	Vector direction;

};


// ====================================================================================================================
// Geometry abstract class
// ====================================================================================================================
class Geometry {
public:
	Geometry(Vector alb, bool mir, bool trans, bool refract) : albedo(alb), mirror(mir), transparency(trans), refraction_index(refract) {}

	virtual bool intersect(const Ray&, Vector&, Vector&, double&) const = 0;

	Vector albedo;
	bool mirror;
	bool transparency;
	double refraction_index;
};


// ====================================================================================================================
// Class Sphere
// ====================================================================================================================
class Sphere : public Geometry {
public:
	explicit Sphere(
		const Vector& C,
		double r,
		const Vector rho=Vector(1, 1, 1),
		bool mirror=false,
		bool transparency=false,
		double refraction=0.,
		bool inverted_normal=false
		) : Geometry:: Geometry(rho, mirror, transparency, refraction) {
			center = C;
			radius = r;
			inverted_normals = inverted_normal;
		}

	bool intersect(const Ray& ray, Vector& P, Vector& N, double& t) const {
		double a = 1;
		double b = dot(ray.direction, ray.origin - center);
		double c = (ray.origin - center).norm2() - sqr(radius);

		double delta = sqr(b) - a * c;

		if (delta < 0) return false;

		double sqrtDelta = sqrt(delta);

		double t0 = - b - sqrtDelta;
		double t1 = - b + sqrtDelta;

		if (t1 < 0) return false;

		if (t0 < 0) t = t1;
		else t = t0;


		P = ray.origin + t * ray.direction;
		if (inverted_normals) N = center - P;
		else N = P - center;
		N.normalize();

		return true;
	}

	Vector center;
	double radius;
	bool inverted_normals;
};


// ====================================================================================================================
// .obj mesh loading
// ====================================================================================================================
class TriangleIndices {
public:
	TriangleIndices(
		int vtxi = -1,
		int vtxj = -1,
		int vtxk = -1,
		int ni = -1,
		int nj = -1,
		int nk = -1,
		int uvi = -1,
		int uvj = -1,
		int uvk = -1,
		int group = -1,
		bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {}
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


class BoundingBox {
public:
	BoundingBox() {
		double t_min = std::numeric_limits<double>::lowest();
		double t_max = std::numeric_limits<double>::max();
		point_max = Vector(t_min, t_min, t_min);
		point_min = Vector(t_max, t_max, t_max);
	}
	BoundingBox(Vector p_min, Vector p_max) {
		point_min = p_min;
		point_max = p_max;
		centroid = 0.5 * (p_min + p_max);
		Vector extent = point_max - point_min;
		area = 2 * (extent[1] * (extent[0] + extent[2]) + extent[2] * extent[0]);
		if (extent[0] > extent[1]) {
			if (extent[0] > extent[2]) where_max_extent = 0;
			else where_max_extent = 2;
		}
		else {
			if (extent[1] > extent[2]) where_max_extent = 1;
			else where_max_extent = 2;
		}
	}

	bool intersect_ray(const Ray& r) const {
		double tx1 = (point_min[0] - r.origin[0]) / r.direction[0];
		double tx2 = (point_max[0] - r.origin[0]) / r.direction[0];
		double txMin = std::min(tx1, tx2);
		double txMax = std::max(tx1, tx2);

		double ty1 = (point_min[1] - r.origin[1]) / r.direction[1];
		double ty2 = (point_max[1] - r.origin[1]) / r.direction[1];
		double tyMin = std::min(ty1, ty2);
		double tyMax = std::max(ty1, ty2);

		double tz1 = (point_min[2] - r.origin[2]) / r.direction[2];
		double tz2 = (point_max[2] - r.origin[2]) / r.direction[2];
		double tzMin = std::min(tz1, tz2);
		double tzMax = std::max(tz1, tz2);

		return std::min(txMax, std::min(tyMax, tzMax)) > std::max(txMin, std::max(tyMin, tzMin));
	}

	void compute() {
		Vector extent = point_max - point_min;
		centroid = 0.5 * (point_min + point_max);
		area = 2 * (extent[1] * (extent[0] + extent[2]) + extent[2] * extent[0]);
		if (extent[0] > extent[1]) {
			if (extent[0] > extent[2]) where_max_extent = 0;
			else where_max_extent = 2;
		}
		else {
			if (extent[1] > extent[2]) where_max_extent = 1;
			else where_max_extent = 2;
		}
	}


	void merge(const BoundingBox &a) {
		Vector new_p_max, new_p_min;
		new_p_max = Vector(std::max(a.point_max[0], point_max[0]),
						   std::max(a.point_max[1], point_max[1]),
						   std::max(a.point_max[2], point_max[2]));
		new_p_min = Vector(std::min(a.point_min[0], point_min[0]),
						   std::min(a.point_min[1], point_min[1]),
						   std::min(a.point_min[2], point_min[2]));
		point_max = new_p_max;
		point_min = new_p_min;
	}



	Vector point_min, point_max, centroid;
	double area;
	int where_max_extent;
};


BoundingBox Union(const BoundingBox &a, const BoundingBox &b) {
	BoundingBox new_bbox;
	new_bbox.point_max = Vector(std::max(a.point_max[0], b.point_max[0]),
								std::max(a.point_max[1], b.point_max[1]),
								std::max(a.point_max[2], b.point_max[2]));
	new_bbox.point_min = Vector(std::min(a.point_min[0], b.point_min[0]),
								std::min(a.point_min[1], b.point_min[1]),
								std::min(a.point_min[2], b.point_min[2]));
	return new_bbox;
}

struct BboxTree {
	BboxTree () {
		children[0] = children[1] = nullptr;
	}

	void initLeaf(const BoundingBox& bbox, int fTO, int nT) {
		etiquette.point_max[0] = 1;
		etiquette.point_max = bbox.point_max;
		etiquette.point_min = bbox.point_min;
		splitAxis = -1;
		firstTriangleOffset = fTO;
		nTriangles = nT;
		children[0] = children[1] = nullptr;
	}

	void initInteriorNode(int sA, BboxTree *child_left, BboxTree *child_right) {
		splitAxis = sA;
 		firstTriangleOffset = -1;
		nTriangles = -1;
		etiquette = Union(child_left->etiquette, child_right->etiquette);
		children[0] = child_left;
		children[1] = child_right;
	}

	int splitAxis, firstTriangleOffset, nTriangles;
	BoundingBox etiquette;
	BboxTree *children[2];
};

struct PointerLessBboxTreeNode {
	PointerLessBboxTreeNode(BoundingBox b) {
		bbox = b;
		children[0] = -1;
		children[1] = -1;
	}

	BoundingBox bbox;
	int firstTriangleOffset, nTriangles;
	int children[2];
};


class TriangleMesh : public Geometry {
public:
	~TriangleMesh() {}
	TriangleMesh(Vector alb = Vector(1, 1, 1), bool mir = false, bool trans = false, bool refract = 1) : Geometry::Geometry(alb, mir, trans, refract) {}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	bool intersect(const Ray& ray, Vector& P, Vector& N, double& t) const {

		bool found_any_intersection = false;

		std::pair<bool, std::vector<std::pair<int, int>>> intersction = bvhRayIntersect(ray, 0);

		if (!intersction.first) return false;

		for (int i = 0; i < intersction.second.size(); i++) {
			for (int j = 0; j < intersction.second[i].second; j++) {
				int offset = intersction.second[i].first;

				Vector A = vertices[orderedTriangles[offset + j]->vtxi];
				Vector B = vertices[orderedTriangles[offset + j]->vtxj];
				Vector C = vertices[orderedTriangles[offset + j]->vtxk];

				Vector e1 = B - A;
				Vector e2 = C - A;
				Vector local_N = cross(e1, e2);
				Vector OA = A - ray.origin;
				Vector OAxU = cross(OA, ray.direction);

				double denom = 1 / dot(ray.direction, local_N);

				double local_t = dot(OA, local_N) * denom;
				if (local_t < 0 or local_t > t)  continue;

				double beta = dot(e2, OAxU) * denom;
				if (beta < 0 or beta > 1) continue;

				double gamma = - dot(e1, OAxU) * denom;
				if (gamma < 0 or gamma > 1) continue;

				double alpha = 1 - beta - gamma;
				if (alpha < 0) continue;

				t = local_t;
				P = A + beta * e1 + gamma * e2;
				N = normals[orderedTriangles[offset + j]->ni] +
					normals[orderedTriangles[offset + j]->nj] * beta +
					normals[orderedTriangles[offset + j]->nk] * gamma;
				found_any_intersection = true;
			}

		}
		if (found_any_intersection) N.normalize();
		return found_any_intersection;
	}


	// bool intersect(const Ray& ray, Vector& P, Vector& N, double& t) const {

	// 	bool found_any_intersection = false;

	// 	if (!bbox.intersect_ray(ray)) return false;

	// 	for (int i = 0; i < indices.size(); i++) {

	// 			Vector A = vertices[indices[i].vtxi];
	// 			Vector B = vertices[indices[i].vtxj];
	// 			Vector C = vertices[indices[i].vtxk];

	// 			Vector e1 = B - A;
	// 			Vector e2 = C - A;
	// 			Vector local_N = cross(e1, e2);
	// 			Vector OA = A - ray.origin;
	// 			Vector OAxU = cross(OA, ray.direction);

	// 			double denom = 1 / dot(ray.direction, local_N);

	// 			double local_t = dot(OA, local_N) * denom;
	// 			if (local_t < 0 or local_t > t)  continue;

	// 			double beta = dot(e2, OAxU) * denom;
	// 			if (beta < 0 or beta > 1) continue;

	// 			double gamma = - dot(e1, OAxU) * denom;
	// 			if (gamma < 0 or gamma > 1) continue;

	// 			double alpha = 1 - beta - gamma;
	// 			if (alpha < 0) continue;

	// 			t = local_t;
	// 			P = A + beta * e1 + gamma * e2;
	// 			N = local_N;
	// 			found_any_intersection = true;


	// 	}
	// 	if (found_any_intersection) N.normalize();
	// 	return found_any_intersection;
	// }


	void transform(double scale_factor, Vector translation, Vector rotation = Vector()) {
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * scale_factor + translation;
		}
	}

	void compute_bbox() {
		double min_d = std::numeric_limits<double>::lowest();
		double max_d = std::numeric_limits<double>::max();
		bbox.point_min = Vector(max_d, max_d, max_d);
		bbox.point_max = Vector(min_d, min_d, min_d);

		for (int i = 0; i < vertices.size(); i++) {
			for (int j = 0; j <3 ; j++) {
				if (vertices[i][j] < bbox.point_min[j]) bbox.point_min[j] = vertices[i][j];
				if (vertices[i][j] > bbox.point_max[j]) bbox.point_max[j] = vertices[i][j];
			}
		}

	}


	void compute_BVH() {
		std::vector<std::pair<int, BoundingBox>> triangleInfo(indices.size());

		for (int i = 0; i < indices.size(); i ++) {
			Vector A = vertices[indices[i].vtxi];
			Vector B = vertices[indices[i].vtxj];
			Vector C = vertices[indices[i].vtxk];

			BoundingBox bboxI(
				Vector(
					std::min(A[0], std::min(B[0], C[0])),
					std::min(A[1], std::min(B[1], C[1])),
					std::min(A[2], std::min(B[2], C[2]))
				),
				Vector(
					std::max(A[0], std::max(B[0], C[0])),
					std::max(A[1], std::max(B[1], C[1])),
					std::max(A[2], std::max(B[2], C[2]))
				)
			);
			triangleInfo[i] = std::pair<int, const BoundingBox>(i, bboxI);
		}

		int treeSize = 0;
		BboxTree *root;
		root = recursiveBuildBVH(triangleInfo, 0, triangleInfo.size(), treeSize);
		flattenBVHTree(root);
	}

	BboxTree* recursiveBuildBVH(std::vector<std::pair<int, BoundingBox>>& triangleInfo, int start, int end, int& treeSize) {
		BboxTree *node = new BboxTree;
		treeSize++;

		int ntriangles = end - start;
		if (ntriangles == 1) {
			// Create leaf node
			// std::cout << "Right before crash" << std::endl;
			// std::cout << "Bounding box : " << triangleInfo[start].second.point_min << triangleInfo[start].second.point_max << std::endl;
			node->initLeaf(triangleInfo[start].second, orderedTriangles.size(), ntriangles);

			for (int i = start; i < end; i++) {
				orderedTriangles.push_back(&indices[triangleInfo[i].first]);
			}
			return node;
		}

		BoundingBox global_bbox = triangleInfo[start].second;
		for (int i = start + 1; i < end; i++) {
			global_bbox.merge(triangleInfo[i].second);
		}
		global_bbox.compute();

		int splitDim = global_bbox.where_max_extent;

		int mid = (start + end) / 2;

		// Check if somehow the biggest extent is 0 i.e. several triangles have the same centroid
		if (global_bbox.point_max[splitDim] == global_bbox.point_min[splitDim]){
			int firstTriangleOffset = orderedTriangles.size();
			for (int i = start; i < end; i++) {
				orderedTriangles.push_back(&indices[triangleInfo[i].first]);
			}
			node->initLeaf(global_bbox, firstTriangleOffset, ntriangles);
			return node;
		}

		// Split according to approximate SAH
		if (ntriangles < 4) {
			// If less than 4 triangles simply use equally size subsets to split
			std::nth_element(
				&triangleInfo[start], &triangleInfo[mid], &triangleInfo[end - 1] + 1,
				[splitDim](const std::pair<int, BoundingBox> &a, const std::pair<int, BoundingBox> &b) {
					return (a.second.centroid[splitDim] < b.second.centroid[splitDim]);
				}
			);

			node ->initInteriorNode(
				splitDim,
				recursiveBuildBVH(
					triangleInfo,
					start,
					mid,
					treeSize
				),
				recursiveBuildBVH(
					triangleInfo,
					mid,
					end,
					treeSize
				)
			);
			return node;
		}
		int nbuckets = 12;
		struct Bucket
		{
			int count = 0;
			BoundingBox bbox;
		};
		Bucket buckets[nbuckets];

		// Fill buckets
		for (int i = start; i < end; i++) {
			double offset = (triangleInfo[i].second.centroid[splitDim] - global_bbox.point_min[splitDim]) /
							(global_bbox.point_max[splitDim] - global_bbox.point_min[splitDim]);

			int b = nbuckets * offset;
			if (b == nbuckets) b -= 1;
			buckets[b].count ++;
			buckets[b].bbox.merge(triangleInfo[i].second);
		}
		for (int n = 0; n < nbuckets; n++) {
			if (buckets[n].count == 0) buckets[n].bbox.area = 0;
			else buckets[n].bbox.compute();
		}

		// Compute approximate SAH cost for each split
		// Travel forward  and backwards to compute the cost of the buckets and its predecessors
		int cumulative_count_forward, cumulative_count_backward;
		double cost_forward[nbuckets - 1], cost_backward[nbuckets - 1];

		BoundingBox cumulative_bbox_forward = buckets[0].bbox, cumulative_bbox_backward = buckets[nbuckets - 1].bbox;

		cumulative_count_forward = buckets[0].count;
		cumulative_count_backward = buckets[nbuckets - 1].count;

		cost_forward[0] = cumulative_count_forward * cumulative_bbox_forward.area;
		cost_backward[nbuckets - 2] = cumulative_count_backward * cumulative_bbox_backward.area;

		for (int n = 1; n < nbuckets - 1; n ++) {
			// Update bounding boxes
			cumulative_bbox_forward.merge(buckets[n].bbox);
			cumulative_bbox_forward.compute();

			cumulative_bbox_backward.merge(buckets[nbuckets - 1 - n].bbox);
			cumulative_bbox_backward.compute();

			cumulative_count_forward += buckets[n].count;
			cumulative_count_backward += buckets[nbuckets - n - 1].count;

			cost_forward[n] = cumulative_count_forward * cumulative_bbox_forward.area;
			cost_backward[nbuckets - 1 - n - 1] = cumulative_count_backward * cumulative_bbox_backward.area;

		}

		if (start == 2058 and end== 2062) throw std::invalid_argument("Oui");
		int where_min_cost = 0;
		double cost = cost_forward[0] + cost_backward[0];
		for (int n = 1; n < nbuckets - 1; n++) {
			if (cost_forward[n] + cost_backward[n] < cost) {
				cost = cost_forward[n] + cost_backward[n];
				where_min_cost = n;
			}
		}

		int ntriangles_before_split = 0;

		for (int n = 0; n < where_min_cost + 1; n ++){
			ntriangles_before_split += buckets[n].count;
		}

		double leafCost = ntriangles;
		if ((cost / global_bbox.area) + 0.125 < leafCost) {
			std::pair<int, BoundingBox> *pmid = std::partition(
				&triangleInfo[start], &triangleInfo[end - 1] + 1,
				[=] (const std::pair<int, BoundingBox> &pi) {
					double offset = (pi.second.centroid[splitDim] - global_bbox.point_min[splitDim]) /
									(global_bbox.point_max[splitDim] - global_bbox.point_min[splitDim]);

					int b = nbuckets * offset;
					if (b == nbuckets) b = nbuckets - 1;
					return b <= where_min_cost;
				}
			);

			mid = pmid - &triangleInfo[0];

			node->initInteriorNode(
				splitDim,
				recursiveBuildBVH(
					triangleInfo,
					start,
					mid,
					treeSize
				),
				recursiveBuildBVH(
					triangleInfo,
					mid,
					end,
					treeSize
				)
			);
		}
		else {
			node->initLeaf(
				global_bbox,
				orderedTriangles.size(),
				ntriangles
			);
			for (int i = start; i < end; i++) {
				orderedTriangles.push_back(&indices[triangleInfo[i].first]);
			}
		}

		return node;
	}


	void flattenBVHTree(BboxTree *node) {
		PointerLessBboxTreeNode new_node(
			node->etiquette
		);

		boundingVolumeHierarchy.push_back(new_node);
		int current_last = boundingVolumeHierarchy.size() - 1;
		if (node->children[0] != nullptr) {
			boundingVolumeHierarchy[current_last].children[0] = boundingVolumeHierarchy.size();
			flattenBVHTree(node->children[0]);

			boundingVolumeHierarchy[current_last].children[1] = boundingVolumeHierarchy.size();
			flattenBVHTree(node->children[1]);
		}
		else {
			boundingVolumeHierarchy[current_last].firstTriangleOffset = node->firstTriangleOffset;
			boundingVolumeHierarchy[current_last].nTriangles = node->nTriangles;
		}
		delete node;
	}

	std::pair<bool, std::vector<std::pair<int, int>>> bvhRayIntersect(const Ray &r, const int current_node) const {
		const PointerLessBboxTreeNode& node = boundingVolumeHierarchy[current_node];
		if (node.children[0] == -1) {
			std::vector<std::pair<int, int>> offsets_and_number_of_triangles;
			offsets_and_number_of_triangles.push_back(std::pair<int, int>(node.firstTriangleOffset, node.nTriangles));
			return std::pair<bool, std::vector<std::pair<int, int>>>(true, offsets_and_number_of_triangles);
		}
		bool found_intersect_below = false;
		std::vector<std::pair<int, int>> result;
		if (boundingVolumeHierarchy[node.children[0]].bbox.intersect_ray(r)) {
			std::pair<bool, std::vector<std::pair<int, int>>> result_left = bvhRayIntersect(r, node.children[0]);
			if (result_left.first) {
				found_intersect_below = true;
				result.reserve(result.size() + distance(result_left.second.begin(), result_left.second.end()));
				result.insert(result.end(), result_left.second.begin(), result_left.second.end());
			}
		}
		if (boundingVolumeHierarchy[node.children[1]].bbox.intersect_ray(r)) {
			std::pair<bool, std::vector<std::pair<int, int>>> result_right = bvhRayIntersect(r, node.children[1]);
			if (result_right.first) {
				found_intersect_below = true;
				result.reserve(result.size() + distance(result_right.second.begin(), result_right.second.end()));
				result.insert(result.end(), result_right.second.begin(), result_right.second.end());
			}
		}
		return std::pair<bool, std::vector<std::pair<int, int>>>(found_intersect_below, result);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	std::vector<PointerLessBboxTreeNode> boundingVolumeHierarchy;
	std::vector<TriangleIndices *> orderedTriangles;

	BoundingBox bbox;

};

// ====================================================================================================================
// Class Scene
// ====================================================================================================================
class Scene {
public:
	explicit Scene(double I, double n = 1): n1(n), I(I) {}

	void add_object(const Geometry& G) {objects.push_back(&G);}

	bool intersect_ray(const Ray& r, Vector& P, Vector& N, int& id, double& t) {

		double local_t;
		Vector local_P, local_N;

		bool has_inter = false;

		for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->intersect(r, local_P, local_N, local_t)) {
				has_inter = true;
				if (local_t <= t) {
					t = local_t;
					P = local_P;
					N = local_N;
					id = i;
				}
			}
		}

		return has_inter;
	}

	Vector getColour(const Ray& ray, int rebond, bool last_bounce_diffuse=false) {

		if (rebond < 0) return Vector(0., 0., 0.); // Finish recursion when potentially overflowing

		Vector P, N;
		int id;
		double t = std::numeric_limits<double>::max();

		if (intersect_ray(ray, P, N, id, t)) {

			Vector colour;

			if (id == 0) {
				if (last_bounce_diffuse) return Vector(0., 0., 0.);
				else return this->I / (4 * sqr(M_PI) * sqr(dynamic_cast<const Sphere*>(this->objects[0])->radius)) * Vector(1., 1., 1.);
			}

			else if (objects[id]->mirror) {
				Vector omega_r = ray.direction - 2 * dot(ray.direction, N) * N;
				Ray rayon_reflechi(P + 10e-3 * N, omega_r);
				return getColour(rayon_reflechi, rebond - 1);
			}

			else if (objects[id]->transparency) {
				bool outside_in = (dot(ray.direction, N) < 0);
				double n2 = objects[id]->refraction_index;
				Vector N_correct;
				double n1_c, n2_c;

				if (outside_in) {
					N_correct = N;
					n1_c = n1;
					n2_c = n2;
				}
				else {
					N_correct = - N;
					n1_c = n2;
					n2_c = n1;
				}
				double cos_theta_i = dot(ray.direction, N_correct);

				if (sqrt(1 - sqr(cos_theta_i)) < (n2_c/n1_c)){
					// Seulement la refraction prise en compte
					Vector omega_t_T = n1_c / n2_c * (ray.direction - cos_theta_i * N_correct);
					Vector omega_t_N = - N_correct * sqrt(1 - sqr(n1_c / n2_c) * (1 - sqr(cos_theta_i)));

					Ray r_refracte(P - 1e-3 * N_correct, omega_t_N + omega_t_T);
					return getColour(r_refracte, rebond - 1);
				}

				else {
					// Reflection totale
					Vector omega_r = ray.direction - 2 * dot(ray.direction, N_correct) * N_correct;
					Ray rayon_reflechi(P + 1e-3 * N_correct, omega_r);
					return getColour(rayon_reflechi, rebond - 1);
				}

			}

			else {
				Vector albed = objects[id]->albedo;
				// Hard direct lighting
				// Vector NPrime, Pprime;
				// int idprime;

				// Ray Raylum(P + 10e-2 * N, vecLum);
				// bool shadows_possible = intersect_ray(Raylum, Pprime, NPrime, idprime);

				// bool shadows_exist = shadows_possible and (P - Pprime).norm2() <d2Lum;

				// if (!shadows_exist) {
				// 	colour = I * albed * std::max(dot(N, vecLum), 0.) / (4 * sqr(M_PI) * d2Lum);
				// }

				// Soft direct lighting
				double lumRadius = dynamic_cast<const Sphere*>(this->objects[0])->radius;
				Vector lumCenter = dynamic_cast<const Sphere*>(this->objects[0])->center;
				Vector dirLum = lumCenter - P;
				dirLum.normalize();

				Vector C_xprime = random_cos(-dirLum);
				Vector xprime = C_xprime * lumRadius + lumCenter;
				double p_xprime = std::max(dot(- dirLum, C_xprime), 1e-8) / (M_PI * sqr(lumRadius));

				Vector x_xprime = xprime - P;
				double d2Lum = x_xprime.norm2();
				x_xprime.normalize();

				Ray Raylum(P + 1e-2 * N, x_xprime);

				Vector NPrime, Pprime;
				int idprime;
				double tprime = 1e15;

				bool shadows_possible = intersect_ray(Raylum, Pprime, NPrime, idprime, tprime);

				bool shadows_exist = shadows_possible and sqr(tprime) < d2Lum * 0.95;

				if (!shadows_exist) {
					colour = I / (4 * sqr(M_PI) * sqr(lumRadius))
					* albed / M_PI * std::max(dot(N, x_xprime), 0.) * std::max(dot(C_xprime, - x_xprime), 0.)
					/ (d2Lum * p_xprime);
				}

				// Compute indirect lighting component
				Vector omegaI = random_cos(N);
				Vector indirectColour = albed * getColour(Ray(P + 1e-3 * N, omegaI), rebond - 1, true);

				colour += indirectColour;

				return colour;
			}
		}
		return Vector();
	}

	std::vector<const Geometry*> objects;

	double I;

	double n1;
};


// ====================================================================================================================
// Main
// ====================================================================================================================
int main() {
	auto t1 = std::chrono::steady_clock::now();

	// paramètre image
    int W = 512;
    int H = 512;

	Sphere S_mirror(Vector(-20, 0, 0), 10, Vector(0.7, 0.5, 0.3), true);
	Sphere S_full_transparent(Vector(0, 0, 0), 10, Vector(0.7, 0.5, 0.3), false, true, 1.4);

	Sphere S_hollow_out(Vector(20, 0, 0.), 10, Vector(0.7, 0.5, 0.3), false, true, 1.25);
	Sphere S_hollow_in(Vector(20, 0, 0.), 9.6, Vector(0.7, 0.5, 0.3), false, true, 1.25, true);

	Sphere S_sol(Vector(0, -1000, 0), 1000 - 10, Vector(0, 0, 1));
	Sphere S_plafond(Vector(0, 1000, 0), 1000 - 60, Vector(1, 0, 0));
	Sphere S_murg(Vector(1000, 0, 0), 1000 - 60, Vector(1, 1, 0));
	Sphere S_murd(Vector(-1000, 0, 0), 1000 - 60, Vector(0, 1, 1));
	Sphere S_fond(Vector(0, 0, -1000), 1000 - 60, Vector(0, 1, 0));
	Sphere S_derriere(Vector(0, 0, 1000), 1000 - 60, Vector(1, 0, 1));

	Vector L(-10, 20, 40);
	double I = 1e10;

	Sphere lumiere(L, 10);

	Scene scene(I);

	scene.add_object(lumiere);

	// scene.add_sphere(S_mirror);
	// scene.add_sphere(S_full_transparent);
	// scene.add_sphere(S_hollow_out);
	// scene.add_sphere(S_hollow_in);
	scene.add_object(S_sol);
	scene.add_object(S_plafond);
	scene.add_object(S_murg);
	scene.add_object(S_murd);
	scene.add_object(S_fond);
	// scene.add_object(S_derriere);

	TriangleMesh mesh;
	mesh.readOBJ("cat.obj");
	mesh.transform(0.6, Vector(0, -10, 0));
	// mesh.compute_bbox();
	mesh.compute_BVH();
	auto t2 = std::chrono::steady_clock::now();

	auto duration_s_bvh = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "BVH building took " << duration_s_bvh << " seconds" << std::endl;

	scene.add_object(mesh);

	Vector camera(0, 0, 55);
	double fov = 60 * M_PI / 180;

	std::vector<unsigned char> image(W*H * 3, 0);

	int maxRays = 64;

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		int thread_num = omp_get_thread_num();

		for (int j = 0; j < W; j++) {
			Vector colour;

			// Launch maxRays ray for anti aliasing
			for (int k = 0; k < maxRays; k++) {

                double r1 = uniform(generator[thread_num]);
                double r2 = uniform(generator[thread_num]);
                double rLog = sqrt(-2 * log(r1));
                double sigma = 0.5;

                // Compute perturbatory element for antialiasing
                double gaussX = rLog * cos(2 * M_PI * r2) * sigma;
                double gaussY = rLog * sin(2 * M_PI * r2) * sigma;

                // Compute the ray that passes through the pixel
                Vector u(j - (double)W * 0.5 + gaussX, -i + (double)H * 0.5 + gaussY, -W / (2*tan(fov/2)));
                u.normalize();
                Ray r(camera, u);

                colour += scene.getColour(r, 20);

			}
			colour = (colour * 1/maxRays);

			// Gamma correction
			image[(i*W + j) * 3 + 0] = std::min(pow(colour[0], 0.45), 255.);   // RED
			image[(i*W + j) * 3 + 1] = std::min(pow(colour[1], 0.45), 255.);  // GREEN
			image[(i*W + j) * 3 + 2] = std::min(pow(colour[2], 0.45), 255.);  // BLUE

		}
	}
	auto t3 = std::chrono::steady_clock::now();

	auto duration_s_render = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
	std::cout << "Rendering took " << duration_s_render << " seconds" << std::endl;
	std::string s = "Test_bvh_" + std::to_string(W) + "_" + std::to_string(H) + "_" + std::to_string(maxRays) + ".png";
	char const *pchar = s.c_str();

	stbi_write_png(pchar, W, H, 3, &image[0], 0);

	return 0;
}