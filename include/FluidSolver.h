#ifndef FLUIDSOLVER_H
#define FLUIDSOLVER_H

#include <ngl/Vec3.h>
#include <vector>

#include "Particle.h"

/// @file FluidSolver.h
/// @brief Position Based Fluids solver class, based on the paper http://mmacklin.com/pbf_sig_preprint.pdf
/// @author Teemu Lindborg
/// @version 1.0
/// @date 18.03.2016 Commenting and tidying up
/// Revision History :
///   Started blocking out 08.02.16
///   Implemented the solver and commented the code 17.03.16
/// @todo Make the code more robust

constexpr float m_pi = 3.14159265359f;

// ---------------------------------------------------------------------------------------
/// @class FluidSolver
/// @brief Solver implementing the algorithm defined in the original PBF paper by M. Macklin & M. Müller
// ---------------------------------------------------------------------------------------
class FluidSolver
{
public:
  // ---------------------------------------------------------------------------------------
  /// @brief FluidSolver Default ctor
  // ---------------------------------------------------------------------------------------
  FluidSolver();

  // ---------------------------------------------------------------------------------------
  /// @brief FluidSolver Default dtor
  // ---------------------------------------------------------------------------------------
  ~FluidSolver() {}

  // ---------------------------------------------------------------------------------------
  /// @brief predictPos Predicts the particle's initial position in the frame and updates
  ///                   its velocity based on the gravity and external forces
  /// @param[io] io_p   Particle that's being updated
  /// @param[in] _t     Time step
  // ---------------------------------------------------------------------------------------
  void predictPos(Particle* &io_p, const float &_t);

  // ---------------------------------------------------------------------------------------
  /// @brief computeLambda        Computes the scaling factor for a particle, used for the position update calculations (formula 11)
  /// @param[io] io_particles     Vector holding all the particles
  /// @param[in] _currentParticle Index of the particle currently being updated
  /// @param[in] _neighbors       Array of neighbor indices for the particle
  /// @param[in] _numNeighbors    Amount of neighbors
  // ---------------------------------------------------------------------------------------
  void computeLambda(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors);

  // ---------------------------------------------------------------------------------------
  /// @brief computeDensity       Computes the density of a particle based on the neighboring particles (formula 2)
  /// @param[io] io_particles     Vector holding all the particles
  /// @param[in] _currentParticle Index of the particle currently being updated
  /// @param[in] _neighbors       Array of neighbor indices for the particle
  /// @param[in] _numNeighbors    Amount of neighbors
  // ---------------------------------------------------------------------------------------
  void computeDensity(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors);

  // ---------------------------------------------------------------------------------------
  /// @brief computeVorticityAndXSPH  Computes and adds xsph viscosity and vorticity confiment to the particles velocity (formulas 16 & 17)
  /// @param[io] io_particles         Vector holding all the particles
  /// @param[in] _currentParticle     Index of the particle currently being updated
  /// @param[in] _neighbors           Array of neighbor indices for the particle
  /// @param[in] _numNeighbors        Amount of neighbors
  /// @param[in] _t                   Time step
  // ---------------------------------------------------------------------------------------
  void computeVorticityAndXSPH(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors, float &_t);

  // ---------------------------------------------------------------------------------------
  /// @brief computeArtificialPressure  Computes artificial pressure correction for particle position update
  /// @param[in] _p                     Vector holding the predicted position of the current particle
  /// @param[in] _n                     Vector holding the predicted position of a neighboring particle
  /// @return                           Correction scalar
  // ---------------------------------------------------------------------------------------
  float computeArtificialPressure(ngl::Vec3 &_p, ngl::Vec3 &_n);

  // ---------------------------------------------------------------------------------------
  /// @brief computeDensityKernel Calculates the weight of a neighboring particle using a poly6 kernel
  /// @param[in] _r               Distance between two particles
  /// @return                     Weight of the neighboring particle
  // ---------------------------------------------------------------------------------------
  float computeDensityKernel(const float &_r);

  // ---------------------------------------------------------------------------------------
  /// @brief computeDensityKernelGradient Calculates the gradient of the density kernel using a spiky kernel
  /// @param[in] _p                       Predicted position of the current particle
  /// @param[in] _n                       Predicted position of a neighboring particle
  /// @return                             Gradient vector
  // ---------------------------------------------------------------------------------------
  ngl::Vec3 computeDensityKernelGradient(ngl::Vec3 &_p, ngl::Vec3 &_n);

  // ---------------------------------------------------------------------------------------
  /// @brief calcPositionUpdate     Calculates a position update for a particle
  /// @param[io] io_particles       Vector holding all the particles
  /// @param[in] _currentParticle   Index of the particle currently being updated
  /// @param[in] _neighbors         Array of neighbor indices for the particle
  /// @param[in] _numNeighbors      Amount of neighbors
  /// @return                       The position update
  // ---------------------------------------------------------------------------------------
  ngl::Vec3 calcPositionUpdate(std::vector<Particle *> &io_particles, const unsigned int &_currentParticle, unsigned int *_neighbors, const unsigned int &_numNeighbors);

private:
  // ---------------------------------------------------------------------------------------
  /// @brief m_n Artificial pressure correction exponent
  // ---------------------------------------------------------------------------------------
  int m_n;

  // ---------------------------------------------------------------------------------------
  /// @brief m_inverseRestDensity 1 / rest density
  // ---------------------------------------------------------------------------------------
  float m_inverseRestDensity;

  // ---------------------------------------------------------------------------------------
  /// @brief m_k Small positive constant for formula 13
  // ---------------------------------------------------------------------------------------
  float m_k;

  // ---------------------------------------------------------------------------------------
  /// @brief m_smoothingLength Smoothing kernel distance threshold
  // ---------------------------------------------------------------------------------------
  float m_smoothingLength;

  // ---------------------------------------------------------------------------------------
  /// @brief m_fixedRadius Small fixed radius inside the smoothing kernel
  // ---------------------------------------------------------------------------------------
  float m_fixedRadius;

  // ---------------------------------------------------------------------------------------
  /// @brief m_epsilon Relaxation parameter
  // ---------------------------------------------------------------------------------------
  float m_epsilon;

  // ---------------------------------------------------------------------------------------
  /// @brief m_xsph_c Scale for the XSPH viscosity
  // ---------------------------------------------------------------------------------------
  float m_xsph_c;

  // ---------------------------------------------------------------------------------------
  /// @brief m_polyKernelConstant Precomputed poly6 kernel constant
  // ---------------------------------------------------------------------------------------
  float m_polyKernelConstant;

  // ---------------------------------------------------------------------------------------
  /// @brief m_spikyKernelConstant Precomputed spiky gradient kernel constant
  // ---------------------------------------------------------------------------------------
  float m_spikyKernelConstant;

  // ---------------------------------------------------------------------------------------
  /// @brief m_gravity Vector holding gravity force
  // ---------------------------------------------------------------------------------------
  ngl::Vec3 m_gravity;

protected:

}; // end of FluidSolver

#endif
