#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include <memory>
#include <ngl/VertexArrayObject.h>
#include <ngl/Vec3.h>

/// @file BoundingBox.h
/// @brief Implementation of a simple BoundingBox used as the boundaries of the simulation
/// @author Teemu Lindborg
/// @version 1.0
/// @date 17.03.2016 Commented
/// Revision History :
///   Initial version 15.02.2016
/// @todo Add the possibility to define each corner point for non axis-aligned boxes.

// ---------------------------------------------------------------------------------------
/// @struct Wall
/// @brief Simple wall structure containing the center point, normal of the wall
// ---------------------------------------------------------------------------------------
typedef struct Wall
{
  ngl::Vec3 centre;
  ngl::Vec3 normal;
  float d;
} Wall;

// ---------------------------------------------------------------------------------------
/// @class BoundingBox
/// @brief Simple BoundingBox class to work as the boundaries of the simulation
// ---------------------------------------------------------------------------------------
class BoundingBox
{
public :
  // ---------------------------------------------------------------------------------------
  /// @brief BoundingBox empty ctor
  // ---------------------------------------------------------------------------------------
  BoundingBox() {}

  // ---------------------------------------------------------------------------------------
  /// @brief BoundingBox default ctor
  /// @param[in] _minx Minimum x-coordinate
  /// @param[in] _maxx Maximum x-coordinate
  /// @param[in] _miny Minimum y-coordinate
  /// @param[in] _maxy Maximum x-coordinate
  /// @param[in] _minz Minimum z-coordinate
  /// @param[in] _maxz Maximum x-coordinate
  // ---------------------------------------------------------------------------------------
  BoundingBox(const ngl::Real &_minx, const ngl::Real &_maxx,
              const ngl::Real &_miny, const ngl::Real &_maxy,
              const ngl::Real &_minz, const ngl::Real &_maxz)
  {
    m_minx = _minx;
    m_maxx = _maxx;
    m_miny = _miny;
    m_maxy = _maxy;
    m_minz = _minz;
    m_maxz = _maxz;
  }

  // ---------------------------------------------------------------------------------------
  /// @brief operator =, overloaded = operator for copying the coordinates of another BoundingBox-type object
  /// @param[in] _rhs BoundingBox to copy
  // ---------------------------------------------------------------------------------------
  void operator =(const BoundingBox &_rhs  ) noexcept
  {
    m_minx = _rhs.m_minx;
    m_maxx = _rhs.m_maxx;
    m_miny = _rhs.m_miny;
    m_maxy = _rhs.m_maxy;
    m_minz = _rhs.m_minz;
    m_maxz = _rhs.m_maxz;
  }

  // ---------------------------------------------------------------------------------------
  /// @brief m_minx Min x-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_minx;

  // ---------------------------------------------------------------------------------------
  /// @brief m_maxx Max x-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_maxx;

  // ---------------------------------------------------------------------------------------
  /// @brief m_miny Min y-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_miny;

  // ---------------------------------------------------------------------------------------
  /// @brief m_maxy Max y-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_maxy;

  // ---------------------------------------------------------------------------------------
  /// @brief m_minz Min z-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_minz;

  // ---------------------------------------------------------------------------------------
  /// @brief m_maxz Max z-coordinate
  // ---------------------------------------------------------------------------------------
  ngl::Real m_maxz;

  // ---------------------------------------------------------------------------------------
  /// @brief m_walls Array containing the 6 walls of the bounding box
  // ---------------------------------------------------------------------------------------
  Wall m_walls[6];

  // ---------------------------------------------------------------------------------------
  /// @brief m_vao Vertex Array Object for drawing the outline of the walls
  // ---------------------------------------------------------------------------------------
  std::unique_ptr<ngl::VertexArrayObject> m_vao;

  // ---------------------------------------------------------------------------------------
  /// @brief buildWalls Method to build the walls based on the given min and max
  ///                   coordinates. Also calculates the normals and center points
  ///                   to be used in the collision detection/response calculations.
  // ---------------------------------------------------------------------------------------
  void buildWalls()
  {
    // Bounding box points
    /*    2_____________6
     *   /|            /|
     * 3/_|__________7/ |
     * |  |          |  |
     * |  |          |  |
     * |  |          |  |
     * |  |          |  |
     * |  |0_________|__|4
     * | /           | /
     * |/____________|/
     * 1              5
     */

    // Define the indices for an indexed vao
    const static GLubyte indices[] = {0, 1, 5, 4, 0, 2, 3, 1, 3, 7, 5, 7, 6, 4, 6, 2};

    // Defining each corner point
    ngl::Vec3 p[8];
    p[0] = ngl::Vec3(m_minx, m_miny, m_minz);
    p[1] = ngl::Vec3(m_minx, m_miny, m_maxz);
    p[2] = ngl::Vec3(m_minx, m_maxy, m_minz);
    p[3] = ngl::Vec3(m_minx, m_maxy, m_maxz);
    p[4] = ngl::Vec3(m_maxx, m_miny, m_minz);
    p[5] = ngl::Vec3(m_maxx, m_miny, m_maxz);
    p[6] = ngl::Vec3(m_maxx, m_maxy, m_minz);
    p[7] = ngl::Vec3(m_maxx, m_maxy, m_maxz);

    // Calculating mid points for wall centre determination
    float halfX = (m_maxx + m_minx)/2.f;
    float halfY = (m_maxy + m_miny)/2.f;
    float halfZ = (m_maxz + m_minz)/2.f;

    // Calculating normals etc for each wall
    m_walls[0].normal = (p[5]-p[1]).cross(p[0]-p[1]);
    m_walls[0].normal.normalize();
    m_walls[0].centre = ngl::Vec3(halfX, m_miny, halfZ);
    m_walls[0].d = -(m_walls[0].normal.m_x * m_walls[0].centre.m_x +
                     m_walls[0].normal.m_y * m_walls[0].centre.m_y +
                     m_walls[0].normal.m_z * m_walls[0].centre.m_z);

    m_walls[1].normal = (p[2]-p[0]).cross(p[1]-p[0]);
    m_walls[1].normal.normalize();
    m_walls[1].centre = ngl::Vec3(m_minx, halfY, halfZ);
    m_walls[1].d = -(m_walls[1].normal.m_x * m_walls[1].centre.m_x +
                     m_walls[1].normal.m_y * m_walls[1].centre.m_y +
                     m_walls[1].normal.m_z * m_walls[1].centre.m_z);

    m_walls[2].normal = (p[3]-p[1]).cross(p[5]-p[1]);
    m_walls[2].normal.normalize();
    m_walls[2].centre = ngl::Vec3(halfX, halfY, m_maxz);
    m_walls[2].d = -(m_walls[2].normal.m_x * m_walls[2].centre.m_x +
                     m_walls[2].normal.m_y * m_walls[2].centre.m_y +
                     m_walls[2].normal.m_z * m_walls[2].centre.m_z);

    m_walls[3].normal = (p[7]-p[5]).cross(p[4]-p[5]);
    m_walls[3].normal.normalize();
    m_walls[3].centre = ngl::Vec3(m_maxx, halfY, halfZ);
    m_walls[3].d = -(m_walls[3].normal.m_x * m_walls[3].centre.m_x +
                     m_walls[3].normal.m_y * m_walls[3].centre.m_y +
                     m_walls[3].normal.m_z * m_walls[3].centre.m_z);

    m_walls[4].normal = (p[6]-p[4]).cross(p[0]-p[4]);
    m_walls[4].normal.normalize();
    m_walls[4].centre = ngl::Vec3(halfX, halfY, m_minz);
    m_walls[4].d = -(m_walls[4].normal.m_x * m_walls[4].centre.m_x +
                     m_walls[4].normal.m_y * m_walls[4].centre.m_y +
                     m_walls[4].normal.m_z * m_walls[4].centre.m_z);

    m_walls[5].normal = (p[7]-p[6]).cross(p[2]-p[6]);
    m_walls[5].normal.normalize();
    m_walls[5].centre = ngl::Vec3(halfX, m_maxy, halfZ);
    m_walls[5].d = -(m_walls[5].normal.m_x * m_walls[5].centre.m_x +
                     m_walls[5].normal.m_y * m_walls[5].centre.m_y +
                     m_walls[5].normal.m_z * m_walls[5].centre.m_z);

    // Setting up the VAO object to draw the outlines of the bounding box
    m_vao.reset( ngl::VertexArrayObject::createVOA(GL_LINE_LOOP) );
    m_vao->bind();
    m_vao->setIndexedData(8*sizeof(ngl::Vec3),
                               p[0].m_x,
                               sizeof(indices),
                               &indices[0],
                               GL_UNSIGNED_BYTE, GL_STATIC_DRAW);

    m_vao->setVertexAttributePointer(0, 3, GL_FLOAT, sizeof(ngl::Vec3), 0);

    m_vao->setNumIndices(sizeof(indices));
    m_vao->unbind();
  }
}; // end of BoundingBox

#endif
