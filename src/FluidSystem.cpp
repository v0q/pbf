#include <iostream>
#include <cmath>
#include <chrono>
//#include <omp.h>
#include "FluidSystem.h"

//----------------------------------------------------------------------------------------------------------------------
FluidSystem::FluidSystem() :
  m_bb(-8.f, 6.f, -10.f, 10.f, -6.5f, 2.0f)
{
  // Initialise some of the member variables
  m_solverIterations = 3;
  m_waves = false;
  m_simulate = false;
}

//----------------------------------------------------------------------------------------------------------------------
FluidSystem::~FluidSystem()
{
  // Clean up the particles after finishing
  std::cout << "Cleaning up " << m_particles.size() << " particles\n";
  for(unsigned int i = 0; i < m_particles.size(); ++i)
    delete m_particles[i];
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSystem::init()
{
  // Spawn particles and add the to a vector
  std::cout << "Building the fluid system\n";
  float scale = 0.24f;
  for (int x = 0; x < 8; x++)
  {
    for (int z = 0; z < 8; z++)
    {
      for (int y = 0; y < 16; y++)
      {
        Particle *p = new Particle();
        p->m_pos = ngl::Vec3(-7.5f, -7.f, -6.f) + scale * ngl::Vec3(x, y, z);
        p->m_colour.set(0.f, 0.62745f, 0.690196f);
        m_particles.push_back(p);
      }
    }
  }

  std::cout << m_particles.size() << " particles spawned\n";

  // Call the grid initialisation function passing it the bounding box, amount of particles
  // and how many neighbors each particle can have (user defined)
  m_nns.init(m_bb, m_particles.size(), 150);

  // Build the walls of the bounding box (VAO, normals etc)
  m_bb.buildWalls();
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSystem::execute()
{
  // Draw the bounding box at the beginning
  m_bb.m_vao->bind();
  m_bb.m_vao->draw();
  m_bb.m_vao->unbind();

  if(m_simulate)
  {
    // If the user wants to "simulate waves", move the bounding box max X wall using a sine function and build the vao, normals etc again
    if(m_waves)
    {
      static float sinWave = 0.f;
      sinWave += 0.035f;
      m_bb.m_maxx = 6.f - fabs(std::sin(sinWave)*5.f);
      m_bb.buildWalls();
    }

    // Time step and inverse timestep used for velocity and position calculations
    float timeStep = 0.016f;
    float invTimeStep = 1.f/timeStep;

    // Parallelising the predicted position and velocity calculations
    // #pragma omp parallel for num_threads(omp_get_max_threads())
    for(unsigned int i = 0; i < m_particles.size(); ++i)
    {
      m_solver.predictPos(m_particles[i], timeStep);
      m_particles[i]->m_posUpdate.set(0.f, 0.f, 0.f);
    }

    // Build the grid and neighbor tables based on the predicted positions
    m_nns.buildTable(m_particles);

    // Iterate the solver
    for(unsigned int iter = 0; iter < m_solverIterations; ++iter)
    {
      // Parallelise the lambda calculation
      // #pragma omp parallel for num_threads(omp_get_max_threads())
      for(unsigned int i = 0; i < m_particles.size(); ++i)
      {
        // Get the neighbors for a particle from the computed table
        std::pair<unsigned int *, unsigned int> neighbors = m_nns.getNeighbors(i);

        // Calculate the density constraint
        m_solver.computeLambda(m_particles, i, neighbors.first, neighbors.second);
      }
      // Parallelise the position update calculation
      // #pragma omp parallel for num_threads(omp_get_max_threads())
      for(unsigned int i = 0; i < m_particles.size(); ++i)
      {
        // Get the neighbors for a particle from the computed table
        std::pair<unsigned int *, unsigned int> neighbors = m_nns.getNeighbors(i);

        // Calculate the position update and handle the environment collisions
        m_particles[i]->m_posUpdate = m_solver.calcPositionUpdate(m_particles, i, neighbors.first, neighbors.second);
        handleEnvCollisions(m_particles[i]);
      }
      // Parallelising
      // Add the position updates to the predicted positions
      // #pragma omp parallel for num_threads(omp_get_max_threads())
      for(unsigned int i = 0; i < m_particles.size(); ++i)
      {
        m_particles[i]->m_predPos += m_particles[i]->m_posUpdate;
      }
    }

    // Parallelising
    // #pragma omp parallel for num_threads(omp_get_max_threads())
    for(unsigned int i = 0; i < m_particles.size(); ++i)
    {
      // Calculate the new velocity for each particle based on the old position and the newly predicted position
      m_particles[i]->m_vel = invTimeStep * (m_particles[i]->m_predPos - m_particles[i]->m_pos);

      // Get the neighbors for a particle from the computed table
      // And compute the vorticity and xsph viscosity
      std::pair<unsigned int *, unsigned int> neighbors = m_nns.getNeighbors(i);
      m_solver.computeVorticityAndXSPH(m_particles, i, neighbors.first, neighbors.second, timeStep);

      // Update the position to be the predicted position
      m_particles[i]->m_pos = m_particles[i]->m_predPos;
    }

    // Clean the grid and the neighbor tables
    m_nns.cleanTable();
  }
}

//----------------------------------------------------------------------------------------------------------------------
void FluidSystem::handleEnvCollisions(Particle *io_p)
{
  ngl::Vec3 newPos, newVel;
  float dist;

  // Loop through the walls
  for(int i = 0; i < 6; ++i)
  {
    // Calculate the distance of the particle from the wall
    dist = io_p->m_predPos.m_x * m_bb.m_walls[i].normal.m_x + io_p->m_predPos.m_y * m_bb.m_walls[i].normal.m_y + io_p->m_predPos.m_z * m_bb.m_walls[i].normal.m_z + m_bb.m_walls[i].d - io_p->m_radius;
    if(dist < 0.0) // Penetrates the wall
    {
      // If the particle penetrated the wall, push it out using the distance of penetration and wall normal
      // and update the predicted position accordingly
      float restcoef = 0.5f;
      newPos = io_p->m_predPos - 2.0*dist*m_bb.m_walls[i].normal;
      newVel = -restcoef*(io_p->m_vel.dot(m_bb.m_walls[i].normal) * m_bb.m_walls[i].normal + (io_p->m_vel - io_p->m_vel.dot(m_bb.m_walls[i].normal) * m_bb.m_walls[i].normal));
      io_p->m_predPos = newPos;
      io_p->m_vel = newVel;
    }
  }
}
