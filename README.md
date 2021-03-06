# Daedalus Playground
Collection of open source OpenGL demos, graphics prototypes and physics sims. The goal was to gather cool looking projects, wrap them around with a tiny and dumb CMake environment and make them run on Windows/Ubuntu. Follow the links to find the authors' original code drop and research. Thank you for releasing your work and source code into the wild!

## Real-time Realistic Ocean Lighting - [Project site](http://evasion.imag.fr/Membres/Eric.Bruneton/)
<p align="center"><img src="Bruneton-Ocean/preview/bruneton_ocean.jpg" width="600" /></p>
<p><b>Abstract:</b> Realistic animation and rendering of the ocean is an important aspect for simulators, movies and video games. By nature, the ocean is a difficult problem for Computer Graphics: it is a dynamic system, it combines wave trains at all scales, ranging from kilometric to millimetric. Worse, the ocean is usually viewed at several distances, from very close to the viewpoint to the horizon, increasing the multi-scale issue, and resulting in aliasing problems. The illumination comes from natural light sources (the Sun and the sky dome), is also dynamic, and often underlines the aliasing issues. In this paper, we present a new algorithm for modelling, animation, illumination and rendering of the ocean, in real-time, at all scales and for all viewing distances. Our algorithm is based on a hierarchical representation, combining geometry, normals and BRDF. For each viewing distance, we compute a simplified version of the geometry, and encode the missing details into the normal and the BRDF, depending on the level of detail required. We then use this hierarchical representation for illumination and rendering. Our algorithm runs in real-time, and produces highly realistic pictures and animations.</p>

## Real-time Animation and Rendering of Ocean Whitecaps - [Project site](https://hal.inria.fr/hal-00967078/)
<p align="center"><img src="Dupuy-Whitecaps/preview/dupuy_whitecaps.jpg" width="600" /></p>
<p><b>Abstract:</b> We present a scalable method to procedurally animate and render vast ocean scenes with whitecaps on the GPU. The whitecap coverage on the ocean surface is determined using a wave deformation criteria which can be pre-filtered linearly. This allows us to take advantage of the fast mip-mapping and texture filtering capabilities of modern hardware and produce plausible and anti-aliased images for scales ranging from centimetric to planetary in real time.</p>

## Multi-resolution screen-space ambient occlusion - [Project site](http://sci.utah.edu/~duong/)
<p align="center"><img src="TDHoang-MRSSAO/preview/tdhoang_mrssao.jpg" width="600" /></p>
<p><b>Abstract:</b> We  present  a  new  screen-space  ambient  occlusion  (SSAO)  algorithm  that  improves  on  the  state-of-the-art  SSAO  methods  in  both  speed  and  image  quality.  Our method computes ambient occlusion (AO) for multiple image resolutions and combines them to obtain the final high-resolution  AO.  It  produces  high-quality  AO  that  includes details from small local occluders to low-frequency occlusions from large  faraway occluders. This approach  allows us to use very small sampling kernels at every resolution, and thereby achieve high performance without resorting to random  sampling,  and  therefore  our  results  do  not  suffer from noise and excessive blur, which are common of SSAO. We use bilateral upsampling to interpolate lower-resolution occlusion values to prevent blockiness and occlusion leaks. Compared to existing SSAO methods, our method produces results closer to ray-traced solutions, while running at comparable or higher frame rates than the fastest SSAO method</p>

## A Chebyshev Semi-Iterative Approach for Accelerating Projective and Position-Based Dynamics (Windows only) - [Project site](http://web.cse.ohio-state.edu/~wang.3602/publications.html)
<p align="center"><img src="HuaminWang-Chebyshev/preview/wang_chebyshev.jpg" width="600" /></p>
<p><b>Abstract:</b> In this paper, we study the use of the Chebyshev semi-iterative approach in projective and position-based dynamics. Although projective dynamics is fundamentally nonlinear, its convergence behavior is similar to that of an iterative method solving a linear system. Because of that, we can estimate the “spectral radius” and use it in the Chebyshev approach to accelerate the convergence by at least one order of magnitude, when the global step is handled by the direct solver, the Jacobi solver, or even the Gauss-Seidel solver. Our experiment shows that the combination of the Chebyshev approach and the direct solver runs fastest on CPU, while the combination of the Chebyshev approach and the Jacobi solver outperforms any other combination on GPU, as it is highly compatible with parallel computing. Our experiment further shows position-based dynamics can be accelerated by the Chebyshev approach as well, although the effect is less obvious for tetrahedral meshes. The whole approach is simple, fast, eﬀective, GPU-friendly, and has a small memory cost.</p>

## Descent Methods for Elastic Body Simulation on the GPU (Windows only) - [Project site](http://web.cse.ohio-state.edu/~wang.3602/publications.html)
<p align="center"><img src="HuaminWang-DescentMethods/preview/wang_descentmethods.jpg" width="600" /></p>
<p><b>Abstract:</b> We show that many existing elastic body simulation approaches can be interpreted as descent methods, under a nonlinear optimization framework derived from implicit timeintegration. The key question is how to ﬁnd an eﬀective descent direction with a low computational cost. Based on this concept, we propose a new gradient descent method using Jacobi preconditioning and Chebyshev acceleration. The convergence rate of this method is comparable to that of LBFGS or nonlinear conjugate gradient. But unlike other methods, it requires no dot product operation, making it suitable for GPU implementation. To further improve its convergence and performance, we develop a series of step length adjustment, initialization, and invertible model conversion techniques, all of which are compatible with GPU acceleration. Our experiment shows that the resulting simulator is simple, fast, scalable, memory-eﬃcient, and robust against very large time steps and deformations. It can correctly simulate the deformation behaviors of many elastic materials, as long as their energy functions are second-order diﬀerentiable and their Hessian matrices can be quickly evaluated. For additional speedups, the method can also serve as a complement to other techniques, such as multi-grid.</p>

## A Scalable Galerkin Multigrid Method for Real-time Simulation of Deformable Objects (Windows only) - [Project site](http://tiantianliu.cn/papers/xian2019multigrid/xian2019multigrid.html)
<p align="center"><img src="Xian-GalerkinMG/preview/xian_galerkinmg.jpg" width="600" /></p>
<p><b>Abstract:</b> We propose a simple yet efficient multigrid scheme to simulate high-resolution deformable objects in their full spaces at interactive frame rates. The point of departure of our method is the Galerkin projection which is simple to construct. However, a naive Galerkin multigrid does not scale well for large and irregular grids because it trades-off matrix sparsity for smaller sized linear systems which eventually stops improving the performance. Given that observation, we design our special projection criterion which is based on skinning space coordinates with piecewise constant weights, to make our Galerkin multigrid method scale for high-resolution meshes without suffering from dense linear solves. The usage of skinning space coordinates enables us to reduce the resolution of grids more aggressively, and our piecewise constant weights further ensure us to always deal with reasonably-sparse linear solves. Our projection matrices also help us to manage multi-level linear systems efficiently. Therefore, our method can be applied to different optimization schemes such as Newton's method and Projective Dynamics, pushing the resolution of a real-time simulation to orders of magnitudes higher. Our final GPU implementation outperforms the other state-of-the-art GPU deformable body simulators, enabling us to simulate large deformable objects with hundred thousands of degrees of freedom in real-time.</p>

## Cuda Fluid and Particles Plugins - Vincent Houze (Windows only) - [Project site](https://vincenthouze.com/portfolio/cuda-fluid-and-particles-plugins/)
<p align="center"><img src="VincentHouze-Fluid/Preview/vincenthouze_fluid.jpg" width="600" /></p>

## Ray Tracing Gems - Volume Path Tracer - [Project site](https://github.com/Apress/ray-tracing-gems/tree/master/Ch_28_Ray_Tracing_Inhomogeneous_Volumes/)
<p align="center"><img src="RayTracingGems-VolumePT/preview/raytracinggems_volumept.jpg" width="600" /></p>

## Bullet Physics Samples - Krishna Kumar
<p align="center"><img src="BulletPhysics-PrimitiveShapes/preview/bullet_primitiveshapes.jpg" width="600" /></p>
<p align="center"><img src="BulletPhysics-Cloth/preview/bullet_cloth.jpg" width="600" /></p>

## Spherical harmonics & Smoke simulation - Miles Macklin - [Project site](https://github.com/mmacklin/sandbox)
<p align="center"><img src="MilesMacklin-Sandbox/preview/mmacklin_sandbox.jpg" width="600" /></p>

## Fluid simulation - Wenzel Jakob - [Project site](http://www.mitsuba-renderer.org/misc.html)
<p align="center"><img src="WenzelJakob-Fluid/preview/wenzeljakob_fluid.jpg" width="600" /></p>

## Fast Simulation of Mass-Spring Systems - [Project site](https://www.cs.utah.edu/~ladislav/liu13fast/liu13fast.html)
<p align="center"><img src="FastMassSpring/preview/fastmassspring.jpg" width="600" /></p>
<p><b>Abstract:</b> We describe a scheme for time integration of mass-spring systems that makes use of a solver based on block coordinate descent. This scheme provides a fast solution for classical linear (Hookean) springs. We express the widely used implicit Euler method as an energy minimization problem and introduce spring directions as auxiliary unknown variables. The system is globally linear in the node positions, and the non-linear terms involving the directions are strictly local. Because the global linear system does not depend on run-time state, the matrix can be pre-factored, allowing for very fast iterations. Our method converges to the same final result as would be obtained by solving the standard form of implicit Euler using Newton's method. Although the asymptotic convergence of Newton's method is faster than ours, the initial ratio of work to error reduction with our method is much faster than Newton's. For real-time visual applications, where speed and stability are more important than precision, we obtain visually acceptable results at a total cost per timestep that is only a fraction of that required for a single Newton iteration. When higher accuracy is required, our algorithm can be used to compute a good starting point for subsequent Newton's iteration.</p>

## Compile and build
```
git clone https://github.com/Woking-34/daedalus-playground.git
cd daedalus-playground
mkdir _build
cd _build
cmake .. -G"Visual Studio 15 2017 Win64"
cmake --build . --target bruneton_ocean --config Release
cmake --build . --target dupuy_whitecaps --config Release
...
```
