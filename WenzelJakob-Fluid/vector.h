#if !defined(__VECTOR_H)
#define __VECTOR_H

class Vector3 {
public:
	Float x, y, z;

	inline Vector3(Float _x = 0.0f, Float _y = 0.0f, Float _z = 0.0f)
		: x(_x), y(_y), z(_z) {
	}
	
	inline Vector3 operator+(const Vector3 &v) const {
		return Vector3(x + v.x, y + v.y, z + v.z);
	}

	inline Vector3 operator-(const Vector3 &v) const {
		return Vector3(x - v.x, y - v.y, z - v.z);
	}

	inline Vector3& operator+=(const Vector3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Vector3& operator-=(const Vector3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Vector3 operator*(Float f) const {
		return Vector3(x*f, y*f, z*f);
	}

	inline Vector3 &operator*=(Float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline Vector3 operator-() const {
		return Vector3(-x, -y, -z);
	}

	inline Vector3 operator/(Float f) const {
		Float r = 1.0f / f;
		return Vector3(x * r, y * r, z * r);
	}

	inline Vector3 &operator/=(Float f) {
		Float r = 1.0f / f;
		x *= r;
		y *= r;
		z *= r;
		return *this;
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline Float lengthSquared() const {
		return x*x + y*y + z*z;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}

	inline bool operator==(const Vector3 &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Vector3 &v) const {
		return !operator==(v);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

class Vector2 {
public:
	Float x, y;

	inline Vector2(Float _x = 0.0f, Float _y = 0.0f)
		: x(_x), y(_y) {
	}
	
	inline Vector2 operator+(const Vector2 &v) const {
		return Vector2(x + v.x, y + v.y);
	}

	inline Vector2 operator-(const Vector2 &v) const {
		return Vector2(x - v.x, y - v.y);
	}

	inline Vector2& operator+=(const Vector2 &v) {
		x += v.x; y += v.y;
		return *this;
	}

	inline Vector2& operator-=(const Vector2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	inline Vector2 operator*(Float f) const {
		return Vector2(x*f, y*f);
	}

	inline Vector2 &operator*=(Float f) {
		x *= f;
		y *= f;
		return *this;
	}

	inline Vector2 operator-() const {
		return Vector2(-x, -y);
	}

	inline Vector2 operator/(Float f) const {
		Float r = 1.0f / f;
		return Vector2(x * r, y * r);
	}

	inline Vector2 &operator/=(Float f) {
		Float r = 1.0f / f;
		x *= r; y *= r;
		return *this;
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline Float lengthSquared() const {
		return x*x + y*y;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}

	inline bool operator==(const Vector2 &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const Vector2 &v) const {
		return !operator==(v);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

inline Float dot(const Vector3 &v1, const Vector3 &v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline Float dot(const Vector2 &v1, const Vector2 &v2) {
	return v1.x*v2.x + v1.y*v2.y;
}

inline Vector3 cross(const Vector3 &v1, const Vector3 &v2) {
	return Vector3(
		(v1.y * v2.z) - (v1.z * v2.y), 
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}


typedef Vector3 Point3;
typedef Vector2 Point2;

#endif /* __VECTOR_H */
