#if !defined(__SPHERE_H)
#define __SPHERE_H

struct Sphere {
	Point3 center;
	Float radius;

	/// Construct an empty sphere
	inline Sphere() {
		radius = 0.0f;
	}

	/// Create a sphere from a given center point and a radius
	inline Sphere(const Point3 &pCenter, Float pRadius)
		: center(pCenter), radius(pRadius) {
	}

	inline bool contains(Point3 p) {
		return (p - center).length() < radius;
	}

	/// Calculate the intersection points with the given ray
	inline bool rayIntersect(const Ray3 &ray, Float &nearHit, Float &farHit) const {
		Vector3 originToCenter = center - ray.o;
		Float distToRayClosest = dot(originToCenter, ray.d);
		Float tmp1 = originToCenter.lengthSquared() - radius*radius;

		if (tmp1 <= 0.0f) {
			/* Inside the sphere */
			nearHit = farHit = 
				std::sqrt(distToRayClosest * distToRayClosest - tmp1)
					+ distToRayClosest;
			return true;
		}

		/* Point3s in different direction */
		if (distToRayClosest < 0.0f)
			return false;

		Float sqrOriginToCenterLength = originToCenter.lengthSquared();
		Float sqrHalfChordDist = radius * radius - sqrOriginToCenterLength
			+ distToRayClosest * distToRayClosest;

		if (sqrHalfChordDist < 0) // Miss
			return false;

		// Hit
		Float hitDistance = std::sqrt(sqrHalfChordDist);
		nearHit = distToRayClosest - hitDistance;
		farHit = distToRayClosest + hitDistance;

		if (nearHit == 0)
			nearHit = farHit;

		return true;
	}

	/// Returns a string representation of the sphere
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "Sphere[center = " << center.toString()
			<< ", radius = " << radius << "]";
		return oss.str();
	}
};

#endif /* __SPHERE_H */
