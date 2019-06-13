#if !defined(__FLUIDSOLVER_H)
#define __FLUIDSOLVER_H

#include "common.h"

const int NDIM = 3;

struct Density {
	Float density;
	Float temp;

	inline Density() : density(0), temp(273) {
	}

	inline Density(Float density, Float temp)
		: density(density), temp(temp) {
	}

	inline Density(Float density, Float temp, Vector3 curl, Float curlMagnitude)
		: density(density), temp(temp) {
	}

	inline Density operator+(const Density &v) const {
		return Density(density + v.density, temp + v.temp);
	}

	inline Density operator*(Float f) const {
		return Density(density*f, temp*f);
	}
};

inline std::ostream& operator<<(std::ostream &os, const Density &d) {
	os << "(" << d.density << ", " << d.temp << ")";
	return os;
}

/**
 * Fluid solver implementation
 *
 * @author Wenzel Jakob
 */
class FluidSolver {
public:
	/**
	 * Creates a new fluid solver instance with the given
	 * grid resolution and voxel size
	 */
	FluidSolver(int gridX, int gridY, int gridZ, Float dx);

	/// Free all associated memory
	~FluidSolver();

	/**
	 * Return an interpolated velocity value for the given position
	 */
	inline Vector3 getVelocity(Point3 p) const {
		/* Re-scale to [0, gridX] x [0,gridY] x [0,gridZ] */
		p.x /= m_dx; p.y /= m_dx; p.z /= m_dx;

		return Vector3(
			getVelocityComponent(p.x, p.y - .5f, p.z - .5f, 0),
			getVelocityComponent(p.x - .5f, p.y, p.z - .5f, 1),
			getVelocityComponent(p.x - .5f, p.y - .5f, p.z, 2)
		);
	}
	
	/**
	 * Return an interpolated velocity value for the given position
	 */
	Vector3 getForce(Point3 p) const {
		/* Re-scale to [0, gridX] x [0,gridY] x [0,gridZ] */
		p.x /= m_dx; p.y /= m_dx; p.z /= m_dx;

		return Vector3(
			getForceComponent(p.x, p.y - .5f, p.z - .5f, 0),
			getForceComponent(p.x - .5f, p.y, p.z - .5f, 1),
			getForceComponent(p.x - .5f, p.y - .5f, p.z, 2)
		);
	}
	
	/** Perform a timestep */
	void step(Float dt);
	
	/** Project onto a divergence-free field */
	void project(Float dt);
	
	/** Add forces (buoyancy etc) */
	void calcForces();

	/** Advect the velocity field */
	void advectVelocityField(Float dt);
	
	/** Advect the density field */
	void advectDensityField(Float dt);

	/** Access methods for the visualization */
	inline int getGridSizeX() const { return m_gridX; }
	inline int getGridSizeY() const { return m_gridY; }
	inline int getGridSizeZ() const { return m_gridZ; }
	inline Float getVoxelSize() const { return m_dx; }
	inline const vector<Density> &getDensity() const { return m_d0; }
	inline const Vector &getCurlMagnitude() const { return m_curlMagnitude; }
	inline const Vector &getDivergence() const { return m_divergence; }
	inline const Vector &getPressure() const { return m_pressure; }
	inline const Vector &getVelocity(int i) const { return m_u0[i]; }
	inline const bool *getSolid() const { return m_solid; }
protected:
	/**
	 * Perform a particle tracing step using the RK2 integrator
	 */
	Point3 traceParticle(const Point3 &p, Float dt) const {
		return clip(p, getVelocity(p + getVelocity(p) * (dt * .5f)) * dt);
	}

	inline Point3 clip(Point3 p, Vector3 v) const {
		Ray3 r(p, v);
		Float nearT, farT;
		if (m_aabb.rayIntersect(r, nearT, farT)) {
			if (nearT < 0)
				nearT = farT;
			if (farT < 0) 
				cout << "Internal error while clipping to AABB!" << endl;
			if (farT > 1)
				return p+v;
			cout << "Clipping to AABB" << endl;

			return p + (v*nearT);
		}
		return p+v;
	}

	/**
	 * Return one of the velocity components by performing tri-linear
	 * interpolation over the staggered grid
	 */
	Float getVelocityComponent(Float x, Float y, Float z, int c) const {
		const int i = (int) x, j = (int) y, k = (int) z;
		const int pos=i+j*(m_gridX+1)+k*m_velSlice;

		/* Return zero for positions outside of the grid */
		if (i < 0 || j < 0 || k < 0 || i >= m_gridX || 
			j >= m_gridY || k >= m_gridZ)
			return 0;

		const Float alpha = x-i,
					beta = y-j,
					gamma = z-k,
					A1 = m_u0[c][pos],
					B1 = m_u0[c][pos+1],
					C1 = m_u0[c][pos+m_gridX+1],
					D1 = m_u0[c][pos+m_gridX+2];
		const Float A2 = m_u0[c][pos+m_velSlice],
					B2 = m_u0[c][pos+m_velSlice+1],
					C2 = m_u0[c][pos+m_velSlice+m_gridX+1],
					D2 = m_u0[c][pos+m_velSlice+m_gridX+2];

		return (1-gamma) * ((1-alpha) * (1-beta) * A1 + alpha * (1-beta) * B1
			 + (1-alpha) * beta * C1 + alpha*beta*D1)
			 + gamma * ((1-alpha) * (1-beta) * A2 + alpha * (1-beta) * B2
			 + (1-alpha) * beta * C2 + alpha*beta*D2);
	}

	/**
	 * Return one of the force components by performing tri-linear
	 * interpolation over the staggered grid
	 */
	Float getForceComponent(Float x, Float y, Float z, int c) const {
		const int i = (int) x, j = (int) y, k = (int) z;
		const int pos=i+j*(m_gridX+1)+k*m_velSlice;

		/* Return zero for positions outside of the grid */
		if (i < 0 || j < 0 || k < 0 || i >= m_gridX || 
			j >= m_gridY || k >= m_gridZ)
			return 0;

		const Float alpha = x-i,
					beta = y-j,
					gamma = z-k,
					A1 = m_F[c][pos],
					B1 = m_F[c][pos+1],
					C1 = m_F[c][pos+m_gridX+1],
					D1 = m_F[c][pos+m_gridX+2];
		const Float A2 = m_F[c][pos+m_velSlice],
					B2 = m_F[c][pos+m_velSlice+1],
					C2 = m_F[c][pos+m_velSlice+m_gridX+1],
					D2 = m_F[c][pos+m_velSlice+m_gridX+2];

		return (1-gamma) * ((1-alpha) * (1-beta) * A1 + alpha * (1-beta) * B1
			 + (1-alpha) * beta * C1 + alpha*beta*D1)
			 + gamma * ((1-alpha) * (1-beta) * A2 + alpha * (1-beta) * B2
			 + (1-alpha) * beta * C2 + alpha*beta*D2);
	}
	
	inline void axpy_prod_fast(const Vector &x, Vector &y) const {
		int i;
		for (i=0; i<m_slice; ++i) {
			Float result = m_ADiag[i] * x[i]
				+ m_APlusX[i] * x[i+1]
				+ m_APlusY[i] * x[i+m_gridX]
				+ m_APlusZ[i] * x[i+m_slice];
			if (i-1 >= 0)
				result += m_APlusX[i-1] * x[i-1];
			if (i-m_gridX >= 0)
				result += m_APlusY[i-m_gridX] * x[i-m_gridX];
			if (i-m_slice >= 0)
				result += m_APlusZ[i-m_slice] * x[i-m_slice];

			y[i] = result;
		}
		for (i=m_numPoints-m_slice; i<m_numPoints; ++i) {
			Float result = m_ADiag[i] * x[i]
				+ m_APlusX[i-1] * x[i-1]
				+ m_APlusY[i-m_gridX] * x[i-m_gridX]
				+ m_APlusZ[i-m_slice] * x[i-m_slice];

			if (i+1 < m_numPoints)
				result += m_APlusX[i] * x[i+1];
			if (i+m_gridX < m_numPoints)
				result += m_APlusY[i] * x[i+m_gridX];
			if (i+m_slice < m_numPoints)
				result += m_APlusZ[i] * x[i+m_slice];

			y[i] = result;
		}
		
		int end = m_numPoints-m_slice;
		#pragma omp parallel
		{
			#pragma omp for schedule(static, 1000)
			for (i=m_slice; i<end; ++i) {
				y[i] = m_ADiag[i] * x[i]
					+ m_APlusX[i] * x[i+1]
					+ m_APlusY[i] * x[i+m_gridX]
					+ m_APlusZ[i] * x[i+m_slice]
					+ m_APlusX[i-1] * x[i-1]
					+ m_APlusY[i-m_gridX] * x[i-m_gridX]
					+ m_APlusZ[i-m_slice] * x[i-m_slice];
			}
		}
	}

	inline void precond_solve(const Vector &b, Vector &x) {
		/* Solve the lower triangular system Lq = b */
		for (int i=0; i<m_numPoints; ++i) {
			if (m_solid[i])
				continue;
			Float temp = b[i];
			if (i > 0)
				temp -= m_APlusX[i-1] * m_precond[i-1] * m_cgQ[i-1];
			if (i > m_gridX)
				temp -= m_APlusY[i-m_gridX] * m_precond[i-m_gridX] * m_cgQ[i-m_gridX];
			if (i > m_slice)
				temp -= m_APlusZ[i-m_slice] * m_precond[i-m_slice] * m_cgQ[i-m_slice];
			m_cgQ[i] = temp * m_precond[i];
		}

		/* Solve the upper triangular system L^Tx = q */
		for (int i=m_numPoints-1; i>=0; --i) {
			if (m_solid[i])
				continue;
			Float temp = m_cgQ[i];
			if (i+1 < m_numPoints)
				temp -= m_APlusX[i] * m_precond[i] * x[i+1];
			if (i+m_gridX < m_numPoints)
				temp -= m_APlusY[i] * m_precond[i] * x[i+m_gridX];
			if (i+m_slice < m_numPoints)
				temp -= m_APlusZ[i] * m_precond[i] * x[i+m_slice];
			x[i] = temp * m_precond[i];
		}
	}
	
	inline Density getDensity(Float x, Float y, Float z) const {
		return getDensityCatmullRom(x, y, z);
	}

	Density getDensityTrilinear(Float x, Float y, Float z) const {
		x = (x / m_dx) - .5f;
		y = (y / m_dx) - .5f;
		z = (z / m_dx) - .5f;

		const int i = (int) x, j = (int) y, k = (int) z, 
				  pos=i+j*m_gridX+k*m_slice;

		if (i < 0 || j < 0 || k < 0 || i > m_gridX-1 || 
			j > m_gridY-1 || k > m_gridZ-1)
			return Density();

		const Float alpha = x-i,
					beta = y-j,
					gamma = z-k;

		const Density
					A1 = m_d0[pos],
					B1 = (i+1<m_gridX) ? m_d0[pos+1] : Density(),
					C1 = (j+1<m_gridY) ? m_d0[pos+m_gridX] : Density(),
					D1 = (i+1<m_gridX && j+1<m_gridY) ? m_d0[pos+m_gridX+1] : Density();

		Density A2, B2, C2, D2;
		if (k + 1 < m_gridZ) {
			A2 = m_d0[pos+m_slice];
			B2 = (i+1<m_gridX) ? m_d0[pos+1+m_slice] : Density();
			C2 = (j+1<m_gridY) ? m_d0[pos+m_gridX+m_slice] : Density();
			D2 = (i+1<m_gridX && j+1<m_gridY) ? m_d0[pos+m_gridX+m_slice+1] : Density();
		}

		return  (A1 * ((1-alpha) * (1-beta))
			  +  B1 * (   alpha  * (1-beta))
			  +  C1 * ((1-alpha) *    beta)
			  +  D1 *     alpha  *    beta) * (1-gamma)
		      + (A2 * ((1-alpha) * (1-beta))
			  +  B2 * (   alpha  * (1-beta))
			  +  C2 * ((1-alpha) *    beta)
			  +  D2 *     alpha  *    beta) * gamma;
	}

	inline Density getDensityCatmullRom(Float x, Float y, Float z) const {
		x = (x / m_dx) - .5f;
		y = (y / m_dx) - .5f;
		z = (z / m_dx) - .5f;

		const int i = (int) x, j = (int) y, k = (int) z;

		if (i < 0 || j < 0 || k < 0 || i > m_gridX-1 || 
			j > m_gridY-1 || k > m_gridZ-1)
			return Density();

		const Float alpha = x-i,
					beta = y-j,
					gamma = z-k;

		const Density
			A = (z-1>= 0) ?     getDensityCatmullRomY(x, y, z-1, alpha, beta) : Density(),
			B =                 getDensityCatmullRomY(x, y, z,   alpha, beta),
			C = (z+1<m_gridZ) ? getDensityCatmullRomY(x, y, z+1, alpha, beta) : Density(),
			D = (z+2<m_gridZ) ? getDensityCatmullRomY(x, y, z+2, alpha, beta) : Density();

		const Float gamma2 = gamma*gamma, gamma3 = gamma2*gamma;
		Density d = A*(-0.5f*gamma + gamma2 - 0.5f*gamma3) + 
			   B*(1.0f - gamma2*(5.0f/2.0f) + gamma3*(3.0f/2.0f)) +
			   C*(0.5f*gamma + 2*gamma2 - gamma3*(3.0f/2.0f)) + 
			   D*(-0.5f*gamma2 + 0.5f*gamma3);
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp) 
         || d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - gamma) + C*gamma;
		}
		return d;
	}

	Density getDensityCatmullRomY(int x, int y, int z, Float alpha, Float beta) const {
		const Density
			A = (y-1>=0) ?      getDensityCatmullRomX(x, y-1, z, alpha) : Density(),
			B =                 getDensityCatmullRomX(x, y,   z, alpha),
			C = (y+1<m_gridY) ? getDensityCatmullRomX(x, y+1, z, alpha) : Density(),
			D = (y+2<m_gridY) ? getDensityCatmullRomX(x, y+2, z, alpha) : Density();

		const Float beta2 = beta*beta, beta3 = beta2*beta;
		Density d = A*(-0.5f*beta + beta2 - 0.5f*beta3) + 
			   B*(1.0f - beta2*(5.0f/2.0f) + beta3*(3.0f/2.0f)) +
			   C*(0.5f*beta + 2*beta2 - beta3*(3.0f/2.0f)) + 
			   D*(-0.5f*beta2 + 0.5f*beta3);
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp) 
         || d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - beta) + C*beta;
		}
		return d;
	}

	inline Density getDensityCatmullRomX(int x, int y, int z, Float t) const {
		const int pos=x+y*m_gridX+z*m_slice;
		const Density A = (x-1 >= 0) ? m_d0[pos-1] : Density(),
					  B = m_d0[pos], 
                      C = (x+1<m_gridX) ? m_d0[pos+1] : Density(),
                      D = (x+2<m_gridX) ? m_d0[pos+2] : Density();
		const Float t2 = t*t, t3 = t2*t;
		Density d = A*(-0.5f*t + t2 - 0.5f*t3) + 
			   B*(1.0f - t2*(5.0f/2.0f) + t3*(3.0f/2.0f)) +
			   C*(0.5f*t + 2*t2 - t3*(3.0f/2.0f)) + 
			   D*(-0.5f*t2 + 0.5f*t3);

		/* Switch to trilinear interpolation in the case of an overshoot */
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp) 
         || d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - t) + C*t;
		}
		return d;
	}

	Float CG(const Vector &b, Vector &x, int max_its, Float tol);
	Float PCG(const Vector &b, Vector &x, int max_its, Float tol);
	void test();
protected:
	Sphere m_sphere;

	/* Smoke densities */
	vector<Density> m_d0, m_d1;

	/* Vorticity confinement related */
	vector<Vector3> m_curl, m_vcForce;
	Vector m_curlMagnitude;

	/* Velocities */
	Vector m_u0[NDIM], m_u1[NDIM],  m_F[NDIM];

	/* Vectors for the Hodge decomposition */
	Vector m_divergence, m_pressure;

	/* Temporary vectors for the PCG iteration */
	Vector m_cgP, m_cgW, m_cgZ, m_cgR, m_cgQ;

	/* Pressure matrix */
	Vector m_ADiag, m_APlusX, m_APlusY, m_APlusZ;
	
	/* Incomplete cholesky preconditioner */
	Vector m_precond;

	/* Boolean vector identifying solid cells */
	bool *m_solid;

	/* Grid resolution */
	int m_gridX, m_gridY, m_gridZ;
	int m_numPoints, m_slice, m_velSlice;
	AABB m_aabb;
	Float m_dx, m_invDx;
	Float m_time;
};

#endif /* __FLUIDSOLVER_H */
