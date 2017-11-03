#ifndef NNS_H
#define NNS_H

#include <unordered_map>
#include "BoundingBox.h"
#include "Particle.h"

/// @file NNS.h
/// @brief Nearest neighbor searching class using a uniform grid where the grid size is based on the particle diameter
/// @author Teemu Lindborg
/// @version 1.0
/// @date 17/03/2016 Commenting and tidying up
/// Revision History :
///   Started blocking out 08/02/2016
///   Implemented the grid and nearest neighbor searching ...-17/03/2016
/// @todo Research and implement a more efficient way

// ---------------------------------------------------------------------------------------
/// @class NNS
/// @brief Uniform grid implementation that creates the grid map and builds neighbor tables
///        for each particle, neighbor table parallelised
// ---------------------------------------------------------------------------------------
class NNS
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief NNS Default ctor
  //----------------------------------------------------------------------------------------------------------------------
  NNS() {}

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief ~NNS Default dtor
  //----------------------------------------------------------------------------------------------------------------------
  ~NNS();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief init                 Method to initialise the grid and prepare it for the nns
  /// @param[in] _bb              Bounding box of the simulation
  /// @param[in] _particleCount   Particle count
  /// @param[in] _maxNeighbors    Max neighbors per particle (if exceeded, the particles will be ignored)
  //----------------------------------------------------------------------------------------------------------------------
  void init(const BoundingBox &_bb, const unsigned int &_particleCount, const unsigned int &_maxNeighbors = 60);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief buildTable       Builds the grid and constructs the neighbor table for each particle
  /// @param[in] _particles   Vector containing the particles
  //----------------------------------------------------------------------------------------------------------------------
  void buildTable(const std::vector<Particle *> &_particles);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief cleanTable Cleans the grid and neighbor tables
  //----------------------------------------------------------------------------------------------------------------------
  void cleanTable();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getNeighbors Method to get the neighbors of a particle
  /// @param[in] _pid     Index of the particle in question
  /// @return             Pair containing the pointer to the neighbor table of the particle and the amount of neighbors
  //----------------------------------------------------------------------------------------------------------------------
  std::pair<unsigned int *, unsigned int> getNeighbors(const int &_pid);

private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief buildNeighborTable Builds the neighbor table
  /// @param[in] _particles     Vector containing the particles
  //----------------------------------------------------------------------------------------------------------------------
  void buildNeighborTable(const std::vector<Particle *> &_particles);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getCellX Method to get the cell's X-coordinate
  /// @param[in] _x   x-coordinate of the particle
  /// @return         x-coordinate of the cell based on the particle's postion
  //----------------------------------------------------------------------------------------------------------------------
  int getCellX(const float &_x);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getCellY Method to get the cell's Y-coordinate
  /// @param[in] _y   y-coordinate of the particle
  /// @return         y-coordinate of the cell based on the particle's postion
  //----------------------------------------------------------------------------------------------------------------------
  int getCellY(const float &_y);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getCellZ Method to get the cell's Z-coordinate
  /// @param[in] _z   z-coordinate of the particle
  /// @return         z-coordinate of the cell based on the particle's postion
  //----------------------------------------------------------------------------------------------------------------------
  int getCellZ(const float &_z);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief getCell  Method to get the 1D cell id
  /// @param[in] _x   x-coordinate of the cell
  /// @param[in] _y   y-coordinate of the cell
  /// @param[in] _z   z-coordinate of the cell
  /// @return         Cell id or -1 if it's outside of the bounding box
  //----------------------------------------------------------------------------------------------------------------------
  int getCell(const int &_x, const int &_y, const int &_z);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_particleCount Amount of particles in the system
  //----------------------------------------------------------------------------------------------------------------------
  unsigned int m_particleCount;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_maxNeighbors Maximum amount of neighbors per particle
  //----------------------------------------------------------------------------------------------------------------------
  unsigned int m_maxNeighbors;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_maxParticlesPerCell Maximum amount of particles per cell
  //----------------------------------------------------------------------------------------------------------------------
  unsigned int m_maxParticlesPerCell;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_neighbors Vector containing the neighbors of each particle
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<unsigned int*> m_neighbors;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_numNeighbors Vector containing the number of particles for each particle
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<unsigned int> m_numNeighbors;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_gridCellNumParticles Vector containing the amount of particles in each cell
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<unsigned int> m_gridCellNumParticles;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_grid Map of the cell id to particle indices
  //----------------------------------------------------------------------------------------------------------------------
  std::unordered_map<unsigned int, std::vector<int>> m_grid;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_fixedRadius Fixed radius to search/accept neighbors from
  //----------------------------------------------------------------------------------------------------------------------
  float m_fixedRadius;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_cells Cell-count for each axis
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_cells;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_cellSize Cell-size for each axis
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_cellSize;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_bb Bounding box of the simulation
  //----------------------------------------------------------------------------------------------------------------------
  BoundingBox m_bb;
}; // end of NNS

#endif
