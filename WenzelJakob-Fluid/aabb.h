/** \brief Axis-aligned bounding box data structure
 */
struct AABB {
public:
	Point3 min;
	Point3 max;

	/// Construct an invalid bounding box
	inline AABB() {
		reset();
	}

	/// Create a bounding box from two 3-dimension vectors
	inline AABB(const Point3 &min, const Point3 &max)
		: min(min), max(max) {
	}

	/// Copy-constructor
	inline AABB(const AABB &aabb) 
		: min(aabb.min), max(aabb.max) {
	}

	/// Assignment operator
	inline AABB &operator=(const AABB &aabb) {
		min = aabb.min;
		max = aabb.max;
		return *this;
	}

	/// Mark the bounding box as invalid
	inline void reset() {
		const Float inf = std::numeric_limits<Float>::infinity();
		min = Point3(inf, inf, inf);
		max = Point3(-inf, -inf, -inf);
	}

	/// Calculate the volume of the bounding box
	inline Float getVolume() const {
		Float x = max.x - min.x;
		Float y = max.y - min.y;
		Float z = max.z - min.z;

		return x*x + y*y + z*z;
	}

	/// Clip to another AABB
	inline void clip(const AABB &aabb) {
		min.x = std::max(min.x, aabb.min.x);
		min.y = std::max(min.y, aabb.min.y);
		min.z = std::max(min.z, aabb.min.z);
		max.x = std::min(max.x, aabb.max.x);
		max.y = std::min(max.y, aabb.max.y);
		max.z = std::min(max.z, aabb.max.z);
	}

	/// Return the center point
	inline Point3 getCenter() const {
		return Point3(max + min) / 2.0f;
	}

	/// Calculate the surface area of the bounding box
	inline Float getSurfaceArea() const {
		Vector3 d = max - min;
		return 2.0f * (d.x*d.y + d.x*d.z + d.y*d.z);
	}

	/// Return the axis with the largest corresponding AABB side
	inline int getLargestAxis() const {
		Vector3 d = max - min;
		if (d.x > d.y && d.x > d.z)
			return 0;
		else if (d.y > d.z)
			return 1;
		else
			return 2;
	}

	/// Return whether this bounding box is valid
	inline bool isValid() const {
		return max.x >= min.x && max.y >= min.y && max.z >= min.z;
	}

	/**
	 * Return whether this bounding box covers a non-zero
	 * amount of space
	 */
	inline bool isEmpty() const {
		return max.x <= min.x
			&& max.y <= min.y
			&& max.z <= min.z;
	}


	/// Return the minimum vector of the bounding box
	inline const Point3 &getMinimum() const {
		return min;
	}

	/// Return the maximum vector of the bounding box
	inline const Point3 &getMaximum() const {
		return max;
	}

	/** \brief Calculate the near and far ray-AABB intersection
	 * points (if they exist).
	 */
	inline bool rayIntersect(Ray3 &ray, Float &nearT, Float &farT) const {
		nearT = -std::numeric_limits<Float>::infinity();
		farT  = std::numeric_limits<Float>::infinity();

		/* For each pair of AABB planes */
		for (int i=0; i<3; i++) {
			const Float direction = ray.d[i];
			const Float origin = ray.o[i];
			const Float minVal = min[i], maxVal = max[i];

			if (direction == 0) {
				/* The ray is parallel to the planes */
				if (origin < minVal || origin > maxVal)
					return false;
			} else {
				/* Calculate intersection distances */
				Float t1 = (minVal - origin) / ray.d[i];
				Float t2 = (maxVal - origin) / ray.d[i];

				if (t1 > t2) {
					Float tmp = t1;
					t1 = t2;
					t2 = tmp;
				}

				nearT = std::max(nearT, t1);
				farT = std::min(farT, t2);

				if (nearT > farT)
					return false;
			}
		}
		return true;
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "AABB[min=" << min.toString() << ", max=" << max.toString() << "]";
		return oss.str();
	}
};

