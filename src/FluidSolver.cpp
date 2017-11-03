#include <iostream>
#include <cmath>
//#include <omp.h>
#include "FluidSolver.h"
#include "Particle.h"

//----------------------------------------------------------------------------------------------------------------------
FluidSolver::FluidSolver()
{
  // Initialises all the solver variables used for the calculations
  // Modify these if you want different results (also fix/break the simulation)
  Particle tmp;
  float h;
  m_n = 4;

  m_inverseRestDensity = 1/1000.f;
  m_k = 0.1f;
  m_smoothingLength = h = tmp.m_radius * 5.f;//0.55f;
  m_fixedRadius = 0.3f * m_smoothingLength;
  m_epsilon = 0.0005f;
  m_xsph_c = 0.002f;

  m_polyKernelConstant = 315.f/( 64.f * m_pi * h*h*h*h*h*h*h*h*h);
  m_spikyKernelConstant = -45.f/( m_pi * h*h*h*h*h*h);
  m_gravity.set(0.f, -9.81f, 0.f);
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSolver::predictPos(Particle* &io_p, const float &_t)
{
  // Add the gravity and external forces to the velocity and predict the new position
  // Also reset the external forces
  io_p->m_vel += m_gravity*_t + io_p->m_extForces*_t;
  io_p->m_predPos = io_p->m_pos + _t*io_p->m_vel;
  io_p->m_extForces.set(0.f, 0.f, 0.f);
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSolver::computeLambda(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors)
{
  // Formulas 8 & 11
  float sumGradientLengthSquared = 0;
  float c = 0.f;
  ngl::Vec3 grad_pi_Ci = ngl::Vec3(0, 0, 0);

  // Calculate the density of the particle based on its neighboring particles
  computeDensity(io_particles, _currentParticle, _neighbors, _numNeighbors);

  // Change the colour of the particle based on the density
  float d = io_particles[_currentParticle]->m_density*m_inverseRestDensity;
  ngl::Vec4 color = d * ngl::Vec4(1.f, 1.f - 0.62745f, 1.f - 0.690196f, 1.f);
  color.set( 0.75f - color.m_x, 1.0f - color.m_y, 1.0f - color.m_z, 1.0f);
  io_particles[_currentParticle]->m_colour.set(color);

  // Solve density constraint
  c = io_particles[_currentParticle]->m_density*m_inverseRestDensity - 1.f;
  if(c > 0.f)
  {
    // Parallelise the calculations here
    // #pragma omp parallel for num_threads(omp_get_max_threads())
    for(unsigned int i = 0; i < _numNeighbors; ++i)
    {
      if(_currentParticle == _neighbors[i])
        continue;

      // Implements the formula 8 of the pbf-paper, accumulates the density kernel gradient
      // to be used to determine density constraint
      ngl::Vec3 accumulatedGradient;
      accumulatedGradient = io_particles[_neighbors[i]]->m_mass * computeDensityKernelGradient(io_particles[_currentParticle]->m_predPos, io_particles[_neighbors[i]]->m_predPos);
      accumulatedGradient *= m_inverseRestDensity;

      sumGradientLengthSquared = sumGradientLengthSquared + accumulatedGradient.dot(accumulatedGradient);
      grad_pi_Ci = grad_pi_Ci + accumulatedGradient;
    }

    // u.u = ||u|| * ||u|| * cos 0 = ||u||^2
    sumGradientLengthSquared += grad_pi_Ci.dot(grad_pi_Ci);
    io_particles[_currentParticle]->m_lambda = -c / (sumGradientLengthSquared + m_epsilon);
  }
  else
  {
    io_particles[_currentParticle]->m_lambda = 0.f;
  }
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSolver::computeDensity(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors)
{
  // Initialise density to 0 and calculate the density using
  // the masses and weights of the neighboring particles (formula 2)
  io_particles[_currentParticle]->m_density = 0.f;
  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for(unsigned int i = 0; i < _numNeighbors; ++i)
  {
    if(_currentParticle == _neighbors[i])
      continue;
    io_particles[_currentParticle]->m_density += io_particles[_neighbors[i]]->m_mass * computeDensityKernel((io_particles[_currentParticle]->m_predPos - io_particles[_neighbors[i]]->m_predPos).length());
  }
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSolver::computeVorticityAndXSPH(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors, float &_t)
{
  ngl::Vec3 vorticity, gradVorticity, tmp, xsphV;

  // Implements functions 15, 16 and 17
  // Parallelise the calculations here
  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for(unsigned int i = 0; i < _numNeighbors; ++i)
  {
    // Skip the particle if it's the current particle
    if(_currentParticle == _neighbors[i])
      continue;

    // Calculate the relative velocity and the vector between the two particles
    ngl::Vec3 v_ij = io_particles[_neighbors[i]]->m_vel - io_particles[_currentParticle]->m_vel;
    ngl::Vec3 p_ij = io_particles[_currentParticle]->m_predPos - io_particles[_neighbors[i]]->m_predPos;
    tmp.cross(v_ij, computeDensityKernelGradient(io_particles[_currentParticle]->m_predPos, io_particles[_neighbors[i]]->m_predPos));

    // Accumulate the cross product of the relative velocity and density kernel gradient
    // to the vorticity force
    vorticity += tmp;

    // Add a viscocity force
    if(io_particles[_neighbors[i]]->m_density != 0.f)
      xsphV += v_ij * computeDensityKernel(p_ij.length());
  }
  // Add the accumulated viscosity to the particle's velocity
  io_particles[_currentParticle]->m_vel += m_xsph_c * xsphV;

  // Calculate a gradient vorticity using the spiky kernel and the accumulated vorticity
  float l = vorticity.length();
  if(l != 0.f)
  {
    // #pragma omp parallel for num_threads(omp_get_max_threads())
    for(unsigned int i = 0; i < _numNeighbors; ++i)
    {
      if(_currentParticle == _neighbors[i])
        continue;
      gradVorticity += computeDensityKernelGradient(io_particles[_currentParticle]->m_predPos, io_particles[_neighbors[i]]->m_predPos) * l;
    }
  }

  if(gradVorticity.lengthSquared() != 0.f)
  {
    // If the gradient vorticity "exists", add it to the external forces
    // This is used for higher splashes
    gradVorticity.normalize();
    io_particles[_currentParticle]->m_extForces += gradVorticity.cross(vorticity) * 0.01f;
  }
}

//----------------------------------------------------------------------------------------------------------------------
float FluidSolver::computeArtificialPressure(ngl::Vec3 &_p, ngl::Vec3 &_n)
{
  // Compute the artificial pressure correction factor
  float scorr, tmp;
  scorr = tmp = computeDensityKernel((_p - _n).length()) / computeDensityKernel(m_fixedRadius);
  for(int i = 1; i < m_n; ++i)
    scorr *= tmp;
  return -m_k * scorr;
}

//----------------------------------------------------------------------------------------------------------------------
float FluidSolver::computeDensityKernel(const float &_r)
{
  // Compute the weight of a particle based on the distance using a poly6 kernel
  if( _r < 0.f || _r > m_smoothingLength )
    return 0.f;

  float tmp = m_smoothingLength*m_smoothingLength - _r*_r;
  return m_polyKernelConstant * tmp*tmp*tmp;
}

//----------------------------------------------------------------------------------------------------------------------
ngl::Vec3 FluidSolver::computeDensityKernelGradient(ngl::Vec3 &_p, ngl::Vec3 &_n)
{
  // Compute the kernel gradient and return the vector between the particles multiplied by this weight
  ngl::Vec3 v = (_p - _n);
  float r = v.length();

  if( r > m_smoothingLength )
    return ngl::Vec3(0.f, 0.f, 0.f);

  v.normalize();
  float tmp = m_smoothingLength - r;
  return m_spikyKernelConstant * tmp*tmp * v;
}

//----------------------------------------------------------------------------------------------------------------------
ngl::Vec3 FluidSolver::calcPositionUpdate(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors)
{
  ngl::Vec3 positionUpdate;

  // Looping through the neighboring particles
  // Parallelising this part
  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for(unsigned int i = 0; i < _numNeighbors; ++i)
  {
    if(_currentParticle == _neighbors[i])
      continue;
    // Implements formula 14
    positionUpdate += (io_particles[_currentParticle]->m_lambda + io_particles[_neighbors[i]]->m_lambda + computeArtificialPressure( io_particles[_currentParticle]->m_predPos, io_particles[_neighbors[i]]->m_predPos )) * computeDensityKernelGradient(io_particles[_currentParticle]->m_predPos, io_particles[_neighbors[i]]->m_predPos);
  }

  return m_inverseRestDensity*positionUpdate;
}
