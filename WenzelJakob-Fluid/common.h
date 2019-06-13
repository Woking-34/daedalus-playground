#include <cstdlib>
#include <cstring>
#include <iostream>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <omp.h>

using std::cout;
using std::cerr;
using std::endl;

using namespace boost::numeric::ublas;
typedef float Float;
#define Epsilon 1e-5

#include "vector.h"

typedef coordinate_matrix<Float> SparseMatrix;
typedef vector<Float> Vector;

/** \brief Simple three-dimensional ray class with
   minimum / maximum extent information */
class Ray3 {
public:
	/// Ray origin
	Point3 o;
	/// Minimum range for intersection tests
	mutable Float mint;
	/// Ray direction
	Vector3 d;
	/// Maximum range for intersection tests
	mutable Float maxt;

	/// Construct a new ray
	Ray3() : mint(Epsilon), maxt(std::numeric_limits<Float>::max()) {
	}

	/// Construct a new ray
	Ray3(Point3 _o, Vector3 _d)
		: o(_o), mint(Epsilon),  d(_d), maxt(std::numeric_limits<Float>::max()) {
	}

	/// Return 3d coordinates of a point on the ray
	Point3 operator() (Float t) const { return o + d*t; }

	/// Return a string representation of this ray
	std::string toString() const {
		std::ostringstream oss;
		oss << "Ray3[orig=" << o.toString() << ", dest=" << d.toString() << "]";
		return oss.str();
	}
};

#include "aabb.h"
#include "sphere.h"

