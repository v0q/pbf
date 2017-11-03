#ifndef FLUIDSYSTEM_H
#define FLUIDSYSTEM_H

#include <memory>
#include <vector>
#include "BoundingBox.h"
#include "FluidSolver.h"
#include "NNS.h"
#include "Particle.h"

/// @file FluidSystem.h
/// @brief Fluid system -class encapsulates and plugs together the whole system, handles the creation of the particles
///                            and calls the correct member classes/methods to build a functioning fluid system
/// @author Teemu Lindborg
/// @version 1.0
/// @date 17/03/2016 Commenting and tidying up
/// Revision History :
///   Started blocking out 08/02/16
///   Implemented the system and commented code -17/03/2016
/// @todo Implement a GUI to run the variables in the system

// ---------------------------------------------------------------------------------------
/// @class FluidSystem
/// @brief Class creating the particles and bringing together the solver and grid
// ---------------------------------------------------------------------------------------
class FluidSystem
{
public:
  // ---------------------------------------------------------------------------------------
  /// @brief FluidSystem Default ctor
  // ---------------------------------------------------------------------------------------
  FluidSystem();

  // ---------------------------------------------------------------------------------------
  /// @brief ~FluidSystem Default dtor
  // ---------------------------------------------------------------------------------------
  ~FluidSystem();

  // ---------------------------------------------------------------------------------------
  /// @brief init Initialises the system, creates particles and initialises the NNS class
  // ---------------------------------------------------------------------------------------
  void init();

  // ---------------------------------------------------------------------------------------
  /// @brief execute Executes the simulation loop
  // ---------------------------------------------------------------------------------------
  void execute();

  // ---------------------------------------------------------------------------------------
  /// @brief toggleSimulation Toggles on and off whether to run the simulation or not
  // ---------------------------------------------------------------------------------------
  void toggleSimulation() { m_simulate ^= true; }

  // ---------------------------------------------------------------------------------------
  /// @brief toggleWaves Toggles on and off whether the "wave" machine should run or not
  // ---------------------------------------------------------------------------------------
  void toggleWaves() { m_waves ^= true; }

  // ---------------------------------------------------------------------------------------
  /// @brief getParticles
  /// @return Vector containing pointers to each particle
  // ---------------------------------------------------------------------------------------
  std::vector<Particle *> getParticles() { return m_particles; }

private:
  // ---------------------------------------------------------------------------------------
  /// @brief handleEnvCollisions  Handles the collision of a particle with the bounding box
  /// @param[io] io_p             Particle that's checked for collisions
  // ---------------------------------------------------------------------------------------
  void handleEnvCollisions(Particle *io_p);

  // ---------------------------------------------------------------------------------------
  /// @brief m_solver Solver class
  // ---------------------------------------------------------------------------------------
  FluidSolver m_solver;

  // ---------------------------------------------------------------------------------------
  /// @brief m_sHash Grid & nearest neighbor search class
  // ---------------------------------------------------------------------------------------
  NNS m_nns;

  // ---------------------------------------------------------------------------------------
  /// @brief m_bb Bounding box of the simulation
  // ---------------------------------------------------------------------------------------
  BoundingBox m_bb;

  // ---------------------------------------------------------------------------------------
  /// @brief m_particles Vector holding all the particles of the system
  // ---------------------------------------------------------------------------------------
  std::vector<Particle *> m_particles;

  // ---------------------------------------------------------------------------------------
  /// @brief m_solverIterations Solver iteration count
  // ---------------------------------------------------------------------------------------
  unsigned int m_solverIterations;

  // ---------------------------------------------------------------------------------------
  /// @brief m_simulate Boolean value to determine whether to run the simulation or not
  // ---------------------------------------------------------------------------------------
  bool m_simulate;

  // ---------------------------------------------------------------------------------------
  /// @brief m_waves Boolean value to determine whether to move the bounding box wall to create waves
  // ---------------------------------------------------------------------------------------
  bool m_waves;

protected:

}; // end of FluidSystem

#endif
