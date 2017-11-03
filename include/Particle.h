#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <ngl/Vec4.h>

/// @file Particle.h
/// @brief Simple particle struct
/// @author Teemu Lindborg
/// @version 1.0
/// @date 08/02/16 Initial blocking
/// Revision History :
/// Started blocking out 08/02/16
/// @todo Refining

// ---------------------------------------------------------------------------------------
/// @struct Particle
/// @brief Simple particle struct holding the position, velocity etc
// ---------------------------------------------------------------------------------------
typedef struct Particle
{
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_pos Position of a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_pos;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_predPos Predicted position of a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_predPos;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_posUpdate Calculated position update of a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_posUpdate;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_vel Velocity of a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_vel;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_extForces External forces to a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_extForces;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_colour Colour of a particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec4 m_colour;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_mass Mass of a particle
  //----------------------------------------------------------------------------------------------------------------------
  float m_mass;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_radius Radius of a particle
  //----------------------------------------------------------------------------------------------------------------------
  float m_radius;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_density Density of a particle
  //----------------------------------------------------------------------------------------------------------------------
  float m_density;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief m_lambda Scaling factor of a particle
  //----------------------------------------------------------------------------------------------------------------------
  float m_lambda;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle Default ctor
  /// @param[in] _x Initial x-position
  /// @param[in] _y Initial y-position
  /// @param[in] _z Initial z-position
  /// @param[in] _r Radius of the particle
  //----------------------------------------------------------------------------------------------------------------------
  Particle(const float &_x = 0.f, const float &_y = 0.f, const float &_z = 0.f, const float &_r = 0.125f) :
    m_pos(_x, _y, _z),
    m_radius(_r),
    m_density(0.f),
    m_lambda(0.f)
  {
    float d = _r*2;
    m_mass = d*d*d*1000.f;
    m_colour = ngl::Vec4(1.0f, 1.0f, 1.0f, 1.0f);
  }
} Particle; // end of struct

#endif