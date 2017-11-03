#include <cmath>
#include <iostream>
//#include <omp.h>
#include "NNS.h"

//----------------------------------------------------------------------------------------------------------------------
NNS::~NNS()
{
  // Clean up the neighbor table and grid
  // As most of the member variables are vectors, we don't have to worry about the cleanup
  for(unsigned int i = 0; i < m_particleCount; ++i)
  {
    delete [] m_neighbors[i];
  }
  m_grid.clear();
}

//----------------------------------------------------------------------------------------------------------------------
void NNS::init(const BoundingBox &_bb, const unsigned int &_particleCount, const unsigned int &_maxNeighbors)
{
  // Initialise the grid
  Particle tmp;
  m_fixedRadius = tmp.m_radius * 5.f;

  m_particleCount = _particleCount;
  m_maxNeighbors = _maxNeighbors;
  m_bb = _bb;

  // Set approximate cell sizes to be 3 times the fixed diameter
  // and calculate the exact cell sizes based on the bounding box
  float diam = m_fixedRadius*2;
  float width = m_bb.m_maxx - m_bb.m_minx;
  float height = m_bb.m_maxy - m_bb.m_miny;
  float depth = m_bb.m_maxz - m_bb.m_minz;

  m_cellSize.set( diam/3.f, diam/3.f, diam/3.f );

  // Calculate the exact cell sizes so they'll fill the space
  m_cells.m_x = std::ceil(width/m_cellSize.m_x);
  m_cells.m_y = std::ceil(height/m_cellSize.m_y);
  m_cells.m_z = std::ceil(depth/m_cellSize.m_z);
  m_cellSize.m_x = width/m_cells.m_x;
  m_cellSize.m_y = height/m_cells.m_y;
  m_cellSize.m_z = depth/m_cells.m_z;

  // Get the cell count, estimate maximum particles per cell
  // resize the vectors and allocate memory for the neighbors of each particle
  unsigned int cellCount = (int)(m_cells.m_x * m_cells.m_y * m_cells.m_z);
  m_maxParticlesPerCell = (int)std::ceil((m_cellSize.m_x*m_cellSize.m_y*m_cellSize.m_z) / (tmp.m_radius*tmp.m_radius*tmp.m_radius)) * 2;

  // Resize the vectors and allocate memory for the neighbor tables
  m_gridCellNumParticles.resize(cellCount);
  m_neighbors.resize(m_particleCount);
  m_numNeighbors.resize(m_particleCount);
  for(unsigned int i = 0; i < m_particleCount; ++i)
  {
    m_neighbors[i] = new unsigned int[m_maxNeighbors];
  }

  // Initialise the grid particle counts to 0 and make sure
  // that the grid cell vectors have enough space and are set to -1 for no particles
  for(unsigned int i = 0; i < cellCount; ++i)
  {
    m_gridCellNumParticles[i] = 0;
    for(unsigned int j = 0; j < m_maxParticlesPerCell; ++j)
      m_grid[i] = std::vector<int>(m_maxParticlesPerCell, -1);
  }
}

//----------------------------------------------------------------------------------------------------------------------
void NNS::buildTable(const std::vector<Particle *> &_particles)
{
  // Iterate over the particles and insert them in their respective cells,
  // ignores the particles if the cell id's invalid or the maximum amount of particles
  // per cell has been reached (in theory this should never be the case)
  for(unsigned int i = 0; i < m_particleCount; ++i)
  {
    const unsigned int x = getCellX(_particles[i]->m_pos.m_x);
    const unsigned int y = getCellY(_particles[i]->m_pos.m_y);
    const unsigned int z = getCellZ(_particles[i]->m_pos.m_z);
    const int cell = getCell(x, y, z);
    if(cell != -1)
    {
      if(m_gridCellNumParticles[cell] < m_maxParticlesPerCell)
      {
        m_grid[cell][m_gridCellNumParticles[cell]++] = i;
      }
    }
  }

  // Build the neighbor tables based on the newly built grid
  buildNeighborTable(_particles);
}

//----------------------------------------------------------------------------------------------------------------------
void NNS::cleanTable()
{
  unsigned int cellCount = (int)(m_cells.m_x * m_cells.m_y * m_cells.m_z);

  // Set the grid cell values to -1 and reset the grid cell particle counts to 0
  // Note that we're not cleaning up the neighbor tables as the m_numNeighbors
  // table makes sure that the particle's will never receive "old"/excess data
  // even if the table's not fully cleaned up
  for(unsigned int i = 0; i < cellCount; ++i)
  {
    for(unsigned int j = 0; j < m_gridCellNumParticles[i]; ++j)
      m_grid[i][j] = -1;
    m_gridCellNumParticles[i] = 0;
  }
}

//----------------------------------------------------------------------------------------------------------------------
void NNS::buildNeighborTable(const std::vector<Particle *> &_particles)
{
  // Get neighboring cells up to 2 cells away (each direction)
  int coords[5] = {0, 1, -1, 2, -2};
  // Can be parallelised
#pragma omp parallel default(shared) num_threads(omp_get_max_threads())
  {
    // Loop through the particles
    #pragma omp for schedule(static)
    for(unsigned int a = 0; a < m_particleCount; ++a)
    {
      // Get the current cell coordinates and initialise the amount of neighbors to 0
      const int x = getCellX(_particles[a]->m_pos.m_x);
      const int y = getCellY(_particles[a]->m_pos.m_y);
      const int z = getCellZ(_particles[a]->m_pos.m_z);
      m_numNeighbors[a] = 0;

      // Loop through the current and neighboring cells
      for(int i = 0; i < 5; ++i)
      {
        for(int j = 0; j < 5; ++j)
        {
          for(int k = 0; k < 5; ++k)
          {
            // Get the 1D cell id and make sure it's valid
            const int cell = getCell(x + coords[i], y + coords[j], z + coords[k]);
            if(cell != -1)
            {
              // Get the particle indices in the cell and iterate over them
              std::vector<int> cellData = m_grid[cell];
              for(unsigned int n = 0; n < m_gridCellNumParticles[cell]; ++n)
              {
                // Get the particle index, p should never be -1 due to us keeping track of the amount
                // of particles per cell but just to be sure we check if this is the case and break out
                // Also don't add the current particle to the table
                int p = cellData[n];
                if(p == -1)
                  break;
                if(p == (int)a)
                  continue;
                // Check that the maximum amount of neighbors hasn't been reached
                if(m_numNeighbors[a] < m_maxNeighbors)
                {
                  // Check if the particle is within a fixed radius of the current particle and add it to the list if so
                  if((_particles[a]->m_pos - _particles[p]->m_pos).lengthSquared() < m_fixedRadius*m_fixedRadius)
                  {
                    m_neighbors[a][m_numNeighbors[a]++] = p;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------
std::pair<unsigned int *, unsigned int> NNS::getNeighbors(const int &_pid)
{
  // Return the neighbor table and count as a std::pair based on the particle index
  std::pair<unsigned int *, unsigned int> neighbors(m_neighbors[_pid], m_numNeighbors[_pid]);
  return neighbors;
}

//----------------------------------------------------------------------------------------------------------------------
int NNS::getCell(const int &_x, const int &_y, const int &_z)
{
  // Validate that the cell coordinates are valid
  // and return the 1D cell id
  // returns -1 for invalid cells
  if(_x < 0 || _x >= (int)m_cells.m_x ||
     _y < 0 || _y >= (int)m_cells.m_y ||
     _z < 0 || _z >= (int)m_cells.m_z)
    return -1;
  else
    return _x + _y*(int)m_cells.m_y + _z*(int)m_cells.m_z*(int)m_cells.m_z;
}

//----------------------------------------------------------------------------------------------------------------------
int NNS::getCellX(const float &_x)
{
  // Moving the coordinate to space from 0 till BB-length and calculating which cell the particle's in (x-wise)
  return (int)( (_x - m_bb.m_minx) / m_cellSize.m_x );
}

//----------------------------------------------------------------------------------------------------------------------
int NNS::getCellY(const float &_y)
{
  // Moving the coordinate to space from 0 till BB-length and calculating which cell the particle's in (y-wise)
  return (int)( (_y - m_bb.m_miny) / m_cellSize.m_y );
}

//----------------------------------------------------------------------------------------------------------------------
int NNS::getCellZ(const float &_z)
{
  // Moving the coordinate to space from 0 till BB-length and calculating which cell the particle's in (z-wise)
  return (int)( (_z - m_bb.m_minz) / m_cellSize.m_z );
}
