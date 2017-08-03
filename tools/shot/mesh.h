#ifndef MESH_H
#define MESH_H

#include <vector>
#include <list>
#include <cstddef>

#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif

template <typename T>
struct vec3d
{
  T x,y,z;

  explicit vec3d() : x(0), y(0), z(0) {}
  vec3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  template <typename S>
  vec3d(const vec3d<S>& s) : x(T(s.x)), y(T(s.y)), z(T(s.z)) {}
  
  vec3d<T> operator+(const vec3d<T>& p) const
  {
	 return vec3d<T>(x + p.x, y + p.y, z + p.z);
  }

  vec3d<T> operator-(const vec3d<T>& p) const
  {
	 return vec3d<T>(x - p.x, y - p.y, z - p.z);
  }

  vec3d<T> operator-() const
  {
	 return vec3d<T>(-x, -y, -z);
  }
  
  vec3d<T>& operator+=(const vec3d<T>& p)
  {
	 x += p.x;
	 y += p.y;
	 z += p.z;
	 return *this;
  }
  
  vec3d<T>& operator-=(const vec3d<T>& p)
  {
	 x -= p.x;
	 y -= p.y;
	 z -= p.z;
	 return *this;
  }
  
  template <typename Scalar>
  vec3d<T>& operator*=(const Scalar& s) {
	 x *= s;
	 y *= s;
	 z *= s;
	 return *this;
  }
  
  template <typename Scalar>
  vec3d<T> operator/(const Scalar& v) const
  {
	 return vec3d<T>(x/v, y/v, z/v);
  }
  
  template <typename Scalar>
  vec3d<T>& operator/=(const Scalar& s) {
	 const T i = T(1) / T(s);
	 x *= i;
	 y *= i;
	 z *= i;
	 return *this;
  }
};

template <typename T>
inline T dot_product(const vec3d<T>& v1, const vec3d<T>& v2) 
{
	return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}

template <typename T>
inline vec3d<T> cross_product(const vec3d<T>& v1, const vec3d<T>& v2)
{
  return vec3d<T>(
		(v1.y*v2.z) - (v1.z*v2.y),
		(v1.z*v2.x) - (v1.x*v2.z),
		(v1.x*v2.y) - (v1.y*v2.x)
		);
}

template <typename T>
inline double magnitude(const vec3d<T>& v1)
{
	return sqrt((v1.x*v1.x) + (v1.y*v1.y) + (v1.z*v1.z));
}

template <typename T>
inline double normalize(vec3d<T>& p)
{
  const double n = magnitude(p);
  if (n==0.) {
	 p.x = 0;
	 p.y = 0;
	 p.z = 0;
  }
  else {
	 p.x /= n;
	 p.y /= n;
	 p.z /= n;
  }
  return n;
}

class oriented_point : public vec3d<double> {
   public:
      vec3d<double> n;

      explicit oriented_point() : vec3d<double>() {}
      oriented_point(const vec3d<double>& p) : vec3d<double>(p) {}
      oriented_point(double xx, double yy, double zz) : vec3d<double>(xx,yy,zz) {}
      oriented_point(const oriented_point& p) : vec3d<double>(p), n(p.n) {}
      oriented_point(const vec3d<double>& p, const vec3d<double>& nn) : vec3d<double>(p), n(nn) {}
   };

class mesh_t
{
public:
	
   struct vertex_data {
      vertex_data() : n_tris(0) {}
      oriented_point p;
      size_t n_tris; // no. of triangles this vertex belongs to
      std::list<int> tri_indices; // triangles this vertex belongs to
   };

   struct triangle_data {
      int p0, p1, p2; // indices of the 3 vertices
	  vec3d<double> n;
   };
   
   std::vector<vertex_data> vertices;
   std::vector<triangle_data> triangles;
   
#ifdef USE_FLANN
   flann::Index< flann::L2<double> >* flann_index;
   double* flann_data;
   flann::Matrix<double>* flann_data_mat;
#endif
   
public:

#ifdef USE_FLANN
	mesh_t() : flann_index(0), flann_data(0), flann_data_mat(0) {}
	
	~mesh_t() throw()
	{
	   if (flann_index)
	   {
			delete[] flann_data_mat->ptr();
			delete flann_data_mat;
			delete flann_index;
	   }
	}
#endif
	
	size_t n_points() const { return vertices.size(); }
	
	size_t n_faces() const { return triangles.size(); }

   const oriented_point& get_vertex(int i) const { return vertices.at(i).p; }
   
#ifdef USE_FLANN
   void update_kd_tree()
   {
		if (flann_index)
	   {
			delete[] flann_data_mat->ptr();
			delete flann_data_mat;
			delete flann_index;
	   }
		
		flann_data = new double[3*vertices.size()];
	
		for (size_t i=0; i<vertices.size(); ++i)
		{
			const oriented_point& p = vertices[i].p;
			flann_data[3*i] = p.x;
			flann_data[3*i+1] = p.y;
			flann_data[3*i+2] = p.z;
		}
		
		flann_data_mat = new flann::Matrix<double>(flann_data, vertices.size(), 3);
		
		// NOTE: L2 actually gives the *squared* euclidean distance
		flann_index = new flann::Index< flann::L2<double> >(*flann_data_mat, flann::KDTreeSingleIndexParams());
		flann_index->buildIndex();
   }
#endif
   
   void put_vertices(const std::vector< vec3d<double> >& points)
	{
	   const size_t n_points = points.size();
	   vertices.resize(n_points);

	   // Put the points in the appropriate structure for the kd-tree
	   // and in the vertices collection for future indexed reference

	   for (size_t k=0; k<n_points; ++k)
	   {
		  // indexed vertices
		  vertex_data v;
		  v.p = vec3d<double>(points[k]);
		  v.n_tris = 0;
		  vertices[k] = v;
	   }
	}
	
	void add_triangle(int p0, int p1, int p2)
	{
	   triangle_data td;

	   td.p0 = p0;
	   td.p1 = p1;
	   td.p2 = p2;
	   
	   const vec3d<double>& p3d0 = vertices.at(p0).p;
	   const vec3d<double>& p3d1 = vertices.at(p1).p;
       const vec3d<double>& p3d2 = vertices.at(p2).p;

	   td.n = vec3d<double>( cross_product(p3d1-p3d0, p3d2-p3d0) );
       normalize(td.n);
	   
	   triangles.push_back(td);
	   const int tri_idx = static_cast<int>(triangles.size() - 1);

	   vertices.at(p0).tri_indices.push_back(tri_idx);
	   vertices.at(p0).n_tris++;
	   vertices.at(p1).tri_indices.push_back(tri_idx);
	   vertices.at(p1).n_tris++;
	   vertices.at(p2).tri_indices.push_back(tri_idx);
	   vertices.at(p2).n_tris++;
	}
	
	void calc_normals()
	{
		const size_t np = vertices.size();
		
		for (size_t k = 0; k < np; ++k)
		{
			vertex_data& vert = vertices[k];
			vert.p.n = vec3d<double>(0,0,0);

			std::list<int>::const_iterator it = vert.tri_indices.begin(), it_end = vert.tri_indices.end();
			for (; it!=it_end; ++it)
				vert.p.n += triangles[ *it ].n;

			double mag = magnitude(vert.p.n);
			vert.p.n /= (mag != 0. ? mag : 1.); // some vertices might be outliers

		} // next vertex
	}
	
#ifdef USE_FLANN
	/**
	 * NOTE: * The query point is NOT included in the neighbors list.
	 *       * The function returns SQUARED distances for efficiency reasons.
	 */
	void nearest_neighbors_with_dist(int p, double radius, std::vector<int>& neighs, std::vector<double>& dists)
	{
	   if (!flann_index)
			update_kd_tree();

		const oriented_point& pt = vertices.at(p).p;
		
		double q[3] = {pt.x, pt.y, pt.z};
		flann::Matrix<double> query((double*)&q, 1, 3);
		
		std::vector< std::vector<int> > out_idx;
		std::vector< std::vector<double> > out_dist;
		
		flann_index->radiusSearch(query, out_idx, out_dist, radius*radius, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
		
		if (out_idx.size() != 1 || out_dist.size() != 1)
			std::cout << "[ERROR] mesh_t::nearest_neighbors_with_dist()" << std::endl;
			
		const std::vector<int>& ns = out_idx.front();
		const std::vector<double>& ds = out_dist.front();
		
		if (ns.front() != p)
		{
			std::cout << "[WARNING] mesh_t::nearest_neighbors_with_dist(): The first neighbor should be the query itself (point " << p << ")." << std::endl;
			neighs = ns;
			dists = ds;
		}
		else
		{
			neighs.resize(ns.size() - 1);
			std::copy(ns.begin()+1, ns.end(), neighs.begin());
			
			dists.resize(ns.size() - 1);
			std::copy(ds.begin()+1, ds.end(), dists.begin());
		}
	}
#endif
	
	/**
	 * Get the indices of the triangles neighboring to a given triangle,
	 * i.e., the triangles that share an edge with it (at most 3 if the mesh is manifold).
	 *
	 * WARNING: - This is currently returning repeating indices (FIXME!), please keep this in mind!
	 *          - This code assumes the mesh is manifold.
	 *
	 * @param tri
	 * @param neighs
	 */
	 void get_tri_imm_neighbors(int tri, std::vector<int>& neighs) const
	 {
		const triangle_data& t = triangles.at(tri);
		const vertex_data& p0 = vertices.at(t.p0), p1 = vertices.at(t.p1), p2 = vertices.at(t.p2);
		
		neighs.reserve(3);
		
		// p0,p1 and p0,p2
		for (std::list<int>::const_iterator it=p0.tri_indices.begin(); it!=p0.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p1 || neigh_t.p1 == t.p1 || neigh_t.p2 == t.p1)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p2 || neigh_t.p1 == t.p2 || neigh_t.p2 == t.p2)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
		
		// p1,p0 and p1,p2
		for (std::list<int>::const_iterator it=p1.tri_indices.begin(); it!=p1.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p0 || neigh_t.p1 == t.p0 || neigh_t.p2 == t.p0)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p2 || neigh_t.p1 == t.p2 || neigh_t.p2 == t.p2)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
		
		// p2,p0 and p2,p1
		for (std::list<int>::const_iterator it=p2.tri_indices.begin(); it!=p2.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p0 || neigh_t.p1 == t.p0 || neigh_t.p2 == t.p0)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p1 || neigh_t.p1 == t.p1 || neigh_t.p2 == t.p1)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
	 }
};
   
#endif
