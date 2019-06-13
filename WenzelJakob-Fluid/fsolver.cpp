#include "fsolver.h"
#define CHUNKSIZE 1

FluidSolver::FluidSolver(int gridX, int gridY, int gridZ, Float dx) 
 : m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_dx(dx) {
	int velocityPoints = (gridX+1)*(gridY+1)*(gridZ+1);
	m_numPoints = gridX*gridY*m_gridZ;
	m_slice = gridX*gridY;
	m_velSlice = (gridX+1)*(gridY+1);
	m_time = 0.0f;

	cout << "FluidSolver: Allocating " << (sizeof(Float)*(m_numPoints*15+velocityPoints*NDIM*3)
		+m_numPoints*2*sizeof(Density))/1024 << " KB for a " 
		<< gridX << "x" << gridY << "x" << gridZ << " MAC grid.." << endl;
	m_aabb = AABB(Vector3(0,0,0),
		Vector3(gridX,gridY,gridZ));
	m_invDx = 1.0f / m_dx;

	/* Allocate dense vectors for all quantities, which are 
	   stored on the MAC grid. */
	m_d0.resize(m_numPoints);
	m_d1.resize(m_numPoints);
	m_cgW.resize(m_numPoints);
	m_cgP.resize(m_numPoints);
	m_cgZ.resize(m_numPoints);
	m_cgR.resize(m_numPoints);
	m_cgQ.resize(m_numPoints);
	m_curl.resize(m_numPoints);
	m_vcForce.resize(m_numPoints);
	m_curlMagnitude.resize(m_numPoints);
	m_divergence.resize(m_numPoints);
	m_pressure.resize(m_numPoints);
	m_precond.resize(m_numPoints);
	m_ADiag.resize(m_numPoints);
	m_APlusX.resize(m_numPoints);
	m_APlusY.resize(m_numPoints);
	m_APlusZ.resize(m_numPoints);
	for (int i=0; i<NDIM; ++i) {
		m_F[i].resize(velocityPoints);
		m_u0[i].resize(velocityPoints);
		m_u1[i].resize(velocityPoints);
	}
	m_solid = new bool[m_numPoints];
	m_time = 0.0f;
	memset(m_solid, 0, sizeof(bool)*m_numPoints);

	/* Add an obstacle */
/*	m_sphere = Sphere(Point3((m_gridX/2)*m_dx,m_gridY/2*m_dx,(m_gridZ/2)*m_dx), m_dx * 8);
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				Point3 p(x * m_dx, y*m_dx, z*m_dx);
				if (m_sphere.contains(p))
					m_solid[pos] = true;
			}
		}
	}
*/	
	cout << "FluidSolver: Constructing pressure matrix" << endl;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				bool fluid = !m_solid[pos];
				bool fluidRight = (x != m_gridX-1) && !m_solid[pos+1];
				bool fluidBelow = (y != m_gridY-1) && !m_solid[pos+m_gridX];
				bool fluidBehind = (z != m_gridZ-1) && !m_solid[pos+m_slice];

				if (fluid && fluidRight) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+1] += 1;
					m_APlusX[pos] = -1;
				}

				if (fluid && fluidBelow) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+m_gridX] += 1;
					m_APlusY[pos] = -1;
				}

				if (fluid && fluidBehind) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+m_slice] += 1;
					m_APlusZ[pos] = -1;
				}
			}
		}
	}

	cout << "FluidSolver: Computing the MIC preconditioner" << endl;
	const Float rho = .25f, tau = 0.97f;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (m_solid[pos])
					continue;
				Float termLeft = 0.0f, termAbove = 0.0f, termFront = 0.0f,
				      termLeft2 = 0.0f, termAbove2 = 0.0f, termFront2 = 0.0f;
				if (x > 0 && !m_solid[pos-1]) {
					termLeft  = m_APlusX[pos-1] * m_precond[pos-1];
					termLeft2 = m_APlusX[pos-1] * (m_APlusY[pos-1] + m_APlusZ[pos-1]) * m_precond[pos-1] * m_precond[pos-1];
				}

				if (y > 0 && !m_solid[pos-m_gridX]) {
					termAbove  = m_APlusY[pos-m_gridX] * m_precond[pos-m_gridX];
					termAbove2 = m_APlusY[pos-m_gridX] * (m_APlusX[pos-m_gridX] + m_APlusZ[pos-m_gridX]) * m_precond[pos-m_gridX] * m_precond[pos-m_gridX];
				}

				if (z > 0 && !m_solid[pos-m_slice]) {
					termFront  = m_APlusZ[pos-m_slice] * m_precond[pos-m_slice];
					termFront2 = m_APlusZ[pos-m_slice] * (m_APlusX[pos-m_slice] + m_APlusY[pos-m_slice]) * m_precond[pos-m_slice] * m_precond[pos-m_slice];
				}

				Float e = m_ADiag[pos] - termLeft*termLeft - termAbove*termAbove
				        - termFront*termFront - tau * (termLeft2 + termAbove2 + termFront2);
				if (e < rho * m_ADiag[pos])
					e = m_ADiag[pos];
				m_precond[pos] = 1.0f / std::sqrt(e);
			}
		}
	}
}

void FluidSolver::calcForces() {
	int z;
	std::fill(m_F[0].begin(), m_F[0].end(), 0);
	std::fill(m_F[1].begin(), m_F[1].end(), 0);
	std::fill(m_F[2].begin(), m_F[2].end(), 0);

	/* Calculate the magnitude of the curl */
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (m_solid[pos])
						continue;
					Vector3 velTop     = getVelocity(Point3((x+0.5f)*m_dx, (y+0.0f)*m_dx, (z+0.5f)*m_dx));
					Vector3 velBot     = getVelocity(Point3((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx));
					Vector3 velLeft    = getVelocity(Point3((x+0.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
					Vector3 velRight   = getVelocity(Point3((x+1.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
					Vector3 velFront   = getVelocity(Point3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.0f)*m_dx));
					Vector3 velBack    = getVelocity(Point3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx));

					Vector3 dudx = (velRight - velLeft) / m_dx;
					Vector3 dudy = (velBot - velTop) / m_dx;
					Vector3 dudz = (velBack - velFront) / m_dx;
					Vector3 curl(dudy.z - dudz.y, dudz.x - dudx.z, dudx.y-dudy.x);
					m_curl[pos] = curl;
					m_curlMagnitude[pos] = curl.length();
				}
			}
		}
	}

	/* Add the vorticity confinement force */
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (y==0 || x==0 || z == 0 || x==m_gridX-1 
						|| y==m_gridY-1 || z == m_gridZ-1 || m_solid[pos]) {
						m_curl[pos] = Vector3(0,0,0);
						continue;
					}
					Vector3 N(
						m_curlMagnitude[pos+1]- m_curlMagnitude[pos-1],
						m_curlMagnitude[pos+m_gridX]- m_curlMagnitude[pos-m_gridX],
						m_curlMagnitude[pos+m_slice]- m_curlMagnitude[pos-m_slice]
					);
					Float length = N.length();
					if (length < Epsilon)
						continue;
					m_vcForce[pos] = cross(N / length, m_curl[pos]) * (m_dx * 0.8f);
				}
			}
		}
	}

	/* Add the buoyancy force */
	Float ambient = 273.0f;
/*	int nCells = 0;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (m_solid[pos]) 
					continue;
				ambient += m_d0[pos].temp;
				if (m_d0[pos].temp != 273)
				++nCells;
			}
		}
	}
	ambient /= nCells;
	cout << ambient << endl;
*/
	const Float a = 0.0625f*0.5f, b=0.025f;
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_slice;
			int velIdx = z * m_velSlice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++velIdx, ++pos) {
					if (m_solid[pos])
						continue;

					const Density top   = getDensity((x+.5f)*m_dx, y*m_dx, (z+.5f)*m_dx);

					if (x != 0 && !m_solid[pos-1]) 
						m_F[0][velIdx] += (m_vcForce[pos].x + m_vcForce[pos-1].x)/2.0f;

					if (y!= 0 && !m_solid[pos-m_gridX]) {
						m_F[1][velIdx] -= -a*top.density + b*(top.temp-ambient);
						m_F[1][velIdx] += (m_vcForce[pos].y + m_vcForce[pos-m_gridX].y)/2.0f;
					}

					if (z != 0 && !m_solid[pos-m_slice]) 
						m_F[2][velIdx] += (m_vcForce[pos].z + m_vcForce[pos-m_slice].z)/2.0f;
				}
			}
		}
	}
}

void FluidSolver::step(Float dt) {
/*	Float maxVelocity = 0.0f;
	for (int z=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x) {
				Point3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx);
				Float length = getVelocity(p).length();
				if (length > maxVelocity)
					maxVelocity = length;
			}
		}
	}
	if (maxVelocity < 0.2)
		maxVelocity = 0.2f;

	dt = std::max(dt, 3*m_dx/maxVelocity);
	dt = 0.1f;
*/
	cout << "FluidSolver: step ("  << dt << ") .." << endl;

	project(dt);
	cout << "  + Advecting velocity field .." << endl;
	advectVelocityField(dt);
	cout << "  + Advecting density field .." << endl;
	advectDensityField(dt);
	cout << "  + Calculating forces .." << endl;
	calcForces();

	Float targetTemp = 273+8.4, rateDensity = 10.0f;//, rateTemp = 1000.0f;
/*	int rad = 10;
	for (int z=-rad; z<=rad; ++z) {
		for (int x=-rad; x<=rad; ++x) {
			Float tmp = rad*rad - x*x-z*z;
			if (tmp < 0)
				continue;
			Float dist = std::sqrt(tmp);
			m_d0[(m_gridX/2+x) + (m_gridY-1 - (int)dist) * m_gridX + (m_gridZ/2+z)*m_slice].density = dt*rateDensity;
			m_d0[(m_gridX/2+x) + (m_gridY-1 - (int)dist) * m_gridX + (m_gridZ/2+z)*m_slice].temp = targetTemp;//273+dist*40;
//				(targetTemp - m_d0[(m_gridX/2+x) + (m_gridY-1) * m_gridX + (m_gridZ/2+z)*m_slice].temp) * dist*1;//(1-std::exp(-dt*rateTemp));
		}
	}
*/

	int rad = 4;
	for (int i=-rad; i<=rad; ++i) {
		for (int j=-rad; j<=rad; ++j) {
			int x=m_gridX-1, y = m_gridY-1 - 6+i, z = m_gridZ/2+j;
			Float tmp = rad*rad - i*i-j*j;
			if (tmp < 0)
				continue;
			int dist = -std::max((int) std::sqrt(tmp)-1, 0);
			m_d0[dist + x + y * m_gridX + z * m_slice].density = dt*rateDensity;
			m_d0[dist + x + y * m_gridX + z * m_slice].temp = targetTemp;
			m_F[0][dist + x + y * (m_gridX+1) + z * m_velSlice] = -3.8f;
		}
	}

	int z;
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int velIdx = z * m_velSlice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++velIdx) {
					m_u0[0][velIdx] += dt*m_F[0][velIdx];
					m_u0[1][velIdx] += dt*m_F[1][velIdx];
					m_u0[2][velIdx] += dt*m_F[2][velIdx];
				}
			}
		}
	}
	m_time += dt;
}

void FluidSolver::project(Float dt) {
	/* Calculate the divergence */
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (int z=0; z<m_gridZ; ++z) {
			int velIdx = z * m_velSlice, pos = z * m_slice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++pos, ++velIdx) {
					if (m_solid[pos]) {
						m_divergence[pos] = 0;
						continue;
					}
					m_divergence[pos] = 
						(m_u0[0][velIdx+1]           - m_u0[0][velIdx]
					+ m_u0[1][velIdx+m_gridX+1]   - m_u0[1][velIdx]
					+ m_u0[2][velIdx+m_velSlice]  - m_u0[2][velIdx])
					* m_invDx;
				}
			}
		}
	}

	PCG(m_divergence * (m_dx*m_dx), m_pressure, 100, 1e-4);

	/* Apply the computed gradients */
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (int z=0; z<m_gridZ; ++z) {
			int pos = z * m_slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (m_solid[pos]) 
						continue;
					int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
					if (x < m_gridX-1 && !m_solid[pos + 1]) 
						m_u1[0][velIdx+1] = m_u0[0][velIdx+1] +
							(m_pressure[pos+1] - m_pressure[pos]) * m_invDx;
					if (y < m_gridY-1 && !m_solid[pos + m_gridX]) 
						m_u1[1][velIdx+m_gridX+1] = m_u0[1][velIdx+m_gridX+1] +
							(m_pressure[pos+m_gridX] - m_pressure[pos]) * m_invDx;
					if (z < m_gridZ-1 && !m_solid[pos + m_slice]) 
						m_u1[2][velIdx+m_velSlice] = m_u0[2][velIdx+m_velSlice] +
							(m_pressure[pos+m_slice] - m_pressure[pos]) * m_invDx;
				}
			}
		}
	}

	m_u0[0].swap(m_u1[0]);
	m_u0[1].swap(m_u1[1]);
	m_u0[2].swap(m_u1[2]);
}


void FluidSolver::advectDensityField(Float dt) {
	/* Advect the density field */
	int z;

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int pos = z*m_slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (pos >= m_numPoints) {
						break;
					}
					if (m_solid[pos])
						continue;
					Point3 p((x+.5f)*m_dx, (y+.5f)*m_dx, (z+.5f)*m_dx);
					Point3 p2 = traceParticle(p, -dt);
					m_d1[pos] = getDensity(p2.x, p2.y, p2.z);
				}
			}
		}
	}
	m_d0.swap(m_d1);
}
	
	
void FluidSolver::advectVelocityField(Float dt) {
	int z;

	/* Advect the velocity field */
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, CHUNKSIZE) 
		for (z=0; z<m_gridZ; ++z) {
			int pos = z*m_slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
					if (m_solid[pos]) 
						continue;

					/* Advect X velocities */
					if (x < m_gridX-1 && !m_solid[pos + 1]) {
						Point3 p((x+1.0f)*m_dx, (y+.5f)*m_dx, (z+0.5f)*m_dx);
						Point3 p2 = traceParticle(p, -dt);
						m_u1[0][velIdx+1] = getVelocity(p2).x;
					}

					/* Advect Y velocities */
					if (y < m_gridY-1 && !m_solid[pos + m_gridX]) {
						Point3 p((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx);
						Point3 p2 = traceParticle(p, -dt);
						m_u1[1][velIdx+m_gridX+1] = getVelocity(p2).y;
					}

					/* Advect Z velocities */
					if (z < m_gridZ-1 && !m_solid[pos + m_slice]) {
						Point3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx);
						Point3 p2 = traceParticle(p, -dt);
						m_u1[2][velIdx+m_velSlice] = getVelocity(p2).z;
					}
				}
			}
		}
	}

	m_u0[0].swap(m_u1[0]);
	m_u0[1].swap(m_u1[1]);
	m_u0[2].swap(m_u1[2]);
}

Float FluidSolver::CG(const Vector &b, Vector &x, int max_its, Float tol) {
	Float beta=0, lastRho=0, rho=0, tolSquared = tol*tol;
	int k=0;

	axpy_prod_fast(x, m_cgR);
	m_cgR = b - m_cgR;
	rho = inner_prod(m_cgR, m_cgR);

	while (k<max_its && rho > tolSquared) {
		if (k == 0) {
			noalias(m_cgP) = m_cgR;
		} else {
			beta = rho / lastRho;
			m_cgP = m_cgR + m_cgP * beta;
		}

		axpy_prod_fast(m_cgP, m_cgW);

		const Float alpha = rho / inner_prod(m_cgP, m_cgW);
		noalias(x) += alpha * m_cgP;
		noalias(m_cgR) -= alpha * m_cgW;
		lastRho = rho;
		rho = inner_prod(m_cgR, m_cgR);
		k++;
	}

	cout << "  + CG: Residual after " << k << " iterations : " << std::sqrt(rho) << endl;
	return std::sqrt(rho);
}

Float FluidSolver::PCG(const Vector &b, Vector &x, int max_its, Float tol) {
	Float beta=0, lastRho=0, rho=0, tolSquared = tol*tol;
	int k=0;

	axpy_prod_fast(x, m_cgR);
	m_cgR = b - m_cgR;
	precond_solve(m_cgR, m_cgZ);
	rho = inner_prod(m_cgR, m_cgZ);

	while (k<max_its && rho > tolSquared) {
		if (k == 0) {
			noalias(m_cgP) = m_cgZ;
		} else {
			beta = rho / lastRho;
			m_cgP = m_cgZ + m_cgP * beta;
		}

		axpy_prod_fast(m_cgP, m_cgW);
		const Float alpha = rho / inner_prod(m_cgP, m_cgW);
		noalias(x) += alpha * m_cgP;
		noalias(m_cgR) -= alpha * m_cgW;

		precond_solve(m_cgR, m_cgZ);
		lastRho = rho;
		rho = inner_prod(m_cgR, m_cgZ);
		k++;
	}
	cout << "  + PCG: Residual after " << k << " iterations : " << std::sqrt(rho) << endl;

	return std::sqrt(rho);
}

void FluidSolver::test() {
	Float sum = 0.0f;
	int nCells = 0;
	srand(0);
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z*m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (m_solid[pos]) {
					m_divergence[pos] = 0;
					continue;
				}
				m_divergence[pos] = (Float) rand() / (Float) RAND_MAX; 
				sum += m_divergence[pos];
				++nCells;
			}
		}
	}
	sum /= nCells;
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z*m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (m_solid[pos]) {
					m_divergence[pos] = 0;
					continue;
				}
				m_divergence[pos] -= sum;
			}
		}
	}
	PCG(m_divergence * (m_dx*m_dx), m_pressure, 1000, 1e-6);
//	CG(m_divergence * (m_dx*m_dx), m_pressure, 1000, 1e-6);
}

FluidSolver::~FluidSolver() {
	delete[] m_solid;
}


