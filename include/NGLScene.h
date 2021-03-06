#ifndef NGLSCENE_H__
#define NGLSCENE_H__
#include <ngl/Camera.h>
#include <ngl/Colour.h>
#include <ngl/Light.h>
#include <ngl/Transformation.h>
#include <ngl/Text.h>
#include <QOpenGLWindow>
#include <chrono>
#include "FluidSystem.h"
//----------------------------------------------------------------------------------------------------------------------
/// @file NGLScene.h
/// @brief this class inherits from the Qt OpenGLWindow and allows us to use NGL to draw OpenGL
/// @author Jonathan Macey
/// @version 1.0
/// @date 10/9/13
/// Revision History :
/// This is an initial version used for the new NGL6 / Qt 5 demos
/// @class NGLScene
/// @brief our main glwindow widget for NGL applications all drawing elements are
/// put in this file
//----------------------------------------------------------------------------------------------------------------------

class NGLScene : public QOpenGLWindow
{
  public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor for our NGL drawing class
    /// @param [in] parent the parent window to the class
    //----------------------------------------------------------------------------------------------------------------------
    NGLScene();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor must close down ngl and release OpenGL resources
    //----------------------------------------------------------------------------------------------------------------------
    ~NGLScene();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the initialize class is called once when the window is created and we have a valid GL context
    /// use this to setup any default GL stuff
    //----------------------------------------------------------------------------------------------------------------------
    void initializeGL();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we want to draw the scene
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL();

private:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we resize the window
    //----------------------------------------------------------------------------------------------------------------------
    // Qt 5.5.1 must have this implemented and uses it
    void resizeGL(QResizeEvent *_event);
    // Qt 5.x uses this instead! http://doc.qt.io/qt-5/qopenglwindow.html#resizeGL
    void resizeGL(int _w, int _h);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Qt Event called when a key is pressed
    /// @param [in] _event the Qt event to query for size etc
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *_event);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called every time a mouse is moved
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseMoveEvent (QMouseEvent * _event );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is pressed
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mousePressEvent ( QMouseEvent *_event);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is released
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseReleaseEvent ( QMouseEvent *_event );

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse wheel is moved
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void wheelEvent( QWheelEvent *_event);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief window width
    //----------------------------------------------------------------------------------------------------------------------
    int m_width;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief window height
    //----------------------------------------------------------------------------------------------------------------------
    int m_height;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief timerEvent Method that update()'s called from
    /// @param _event Timer event
    //----------------------------------------------------------------------------------------------------------------------
    void timerEvent(QTimerEvent *_event);

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_rotate Boolean value determining whether the user's rotating the scene or not
    //----------------------------------------------------------------------------------------------------------------------
    bool m_rotate;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_origX Mouse x-position before the move event
    //----------------------------------------------------------------------------------------------------------------------
    int m_origX;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_origY Mouse y-position before the move event
    //----------------------------------------------------------------------------------------------------------------------
    int m_origY;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_spinXFace Amount to rotate around X
    //----------------------------------------------------------------------------------------------------------------------
    float m_spinXFace;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_spinYFace Amount to rotate around Y
    //----------------------------------------------------------------------------------------------------------------------
    float m_spinYFace;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_pbf FluidSystem class containing the simulation
    //----------------------------------------------------------------------------------------------------------------------
    FluidSystem m_pbf;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_cam Simple camera used for view and projection matrices
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Camera m_cam;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_text FPS counter
    //----------------------------------------------------------------------------------------------------------------------
    std::unique_ptr<ngl::Text> m_text;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief m_start/m_end Start and end times of a frame
    //----------------------------------------------------------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
}; // end of NGLScnee

#endif
