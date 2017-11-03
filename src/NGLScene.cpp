#include <QMouseEvent>
#include <QGuiApplication>

#include <ngl/NGLInit.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOPrimitives.h>
#include <iostream>
//#include <omp.h>

#include "NGLScene.h"
#include "FluidSystem.h"

NGLScene::NGLScene()
{
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  setTitle("Position Based Fluids");
  m_rotate = false;
  m_cam.setDefaultCamera();
  m_cam.set(ngl::Vec3(2.f, 8.5f, 16.f), ngl::Vec3(0.f, 0.f, 0.f), ngl::Vec3(0.f, 1.f, 0.f));

//  std::cout << "Setting number of threads to " << omp_get_max_threads() << "\n";
//  omp_set_num_threads(omp_get_max_threads());
}

NGLScene::~NGLScene()
{
  std::cout << "Shutting down NGL, removing VAO's and Shaders\n";
//  m_vao->removeVOA();
}

void NGLScene::resizeGL(QResizeEvent *_event)
{
  m_width = _event->size().width() * devicePixelRatio();
  m_height = _event->size().height() * devicePixelRatio();

  m_cam.setShape(90.f, (m_width / devicePixelRatio())/(m_height / devicePixelRatio()), 0.01f, 20.f);
}
#include "FluidSystem.h"

void NGLScene::resizeGL(int _w , int _h)
{
  m_width = _w * devicePixelRatio();
  m_height = _h * devicePixelRatio();

  m_cam.setShape(90.f, (m_width / devicePixelRatio())/(m_height / devicePixelRatio()), 0.01f, 20.f);
}


void NGLScene::initializeGL()
{
  // we need to initialise the NGL lib which will load all of the OpenGL functions, this must
  // be done once we have a valid GL context but before we call any GL commands. If we dont do
  // this everything will crash
  ngl::NGLInit::instance();
  glClearColor(0.6f, 0.6f, 0.6f, 1.0f);			   // Grey Background
  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);
  // enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);

  ngl::VAOPrimitives *particle = ngl::VAOPrimitives::instance();
  ngl::ShaderLib *shader = ngl::ShaderLib::instance();

  // Creating a VAOPrimitive sphere for particle drawing
  particle->createSphere("particle", 1.f, 32);

  // Create a simple colour shader and set the initial uniforms
  shader->createShaderProgram("SimpleShader");
  shader->attachShader("SimpleVertex", ngl::ShaderType::VERTEX);
  shader->attachShader("SimpleFragment", ngl::ShaderType::FRAGMENT);

  shader->loadShaderSource("SimpleVertex", "shaders/simple.vert");
  shader->loadShaderSource("SimpleFragment", "shaders/simple.frag");

  shader->compileShader("SimpleVertex");
  shader->compileShader("SimpleFragment");

  shader->attachShaderToProgram("SimpleShader", "SimpleVertex");
  shader->attachShaderToProgram("SimpleShader", "SimpleFragment");

  shader->linkProgramObject("SimpleShader");
  shader->use("SimpleShader");

  shader->autoRegisterUniforms("SimpleShader");
  shader->printProperties();

  shader->setRegisteredUniform("u_Light.Position", ngl::Vec4(1.f, 5.5f, 0.5f, 1.f));
  shader->setRegisteredUniform("u_Light.La", ngl::Vec3(0.f, 0.f, 0.f));
  shader->setRegisteredUniform("u_Light.Ld", ngl::Vec3(1.f, 1.f, 1.f));
  shader->setRegisteredUniform("u_Light.Ls", ngl::Vec3(0.1f, 0.1f, 0.1f));

  shader->setRegisteredUniform("u_BackLight.Position", ngl::Vec4(3.f, -1.5f, -25.f, 1.f));
  shader->setRegisteredUniform("u_BackLight.La", ngl::Vec3(0.f, 0.f, 0.f));
  shader->setRegisteredUniform("u_BackLight.Ld", ngl::Vec3(1.f, 1.f, 1.f));
  shader->setRegisteredUniform("u_BackLight.Ls", ngl::Vec3(0.1f, 0.1f, 0.1f));

  shader->setRegisteredUniform("u_Material.Ka", ngl::Vec3(0.2f, 0.2f, 0.2f));
  shader->setRegisteredUniform("u_Material.Kd", ngl::Vec3(1.f, 1.f, 1.f));
  shader->setRegisteredUniform("u_Material.Ks", ngl::Vec3(1.f, 1.f, 1.f));
  shader->setRegisteredUniform("u_Material.Shininess", 2.f);

  // Initialise the fluid system
  m_pbf.init();
  m_text.reset(new ngl::Text(QFont("Arial",14)));
  m_text->setScreenSize(width(),height());

  startTimer(10);
}

void NGLScene::timerEvent(QTimerEvent *_event)
{
  update();
}

void NGLScene::paintGL()
{
  // clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_width,m_height);

  // Get the frame time and set the frame start clock again
  std::chrono::duration<float> elapsed_time = m_end - m_start;
  m_start = std::chrono::system_clock::now();
  std::cout << "Frame time: " << elapsed_time.count() << "s\n";

  // Render out the fps and particle count
  m_text->setColour(1,1,0);
  QString text = QString("%1 fps").arg(1/elapsed_time.count());
  m_text->renderText(10,20,text);
  text=QString("Num particles = %1").arg(m_pbf.getParticles().size());
  m_text->renderText(10,40,text);

  ngl::VAOPrimitives *particle = ngl::VAOPrimitives::instance();
  ngl::ShaderLib *shader = ngl::ShaderLib::instance();

  shader->use("SimpleShader");

  ngl::Mat4 modelMatrix;
  ngl::Mat4 mouseGlobalTX;
  ngl::Mat4 rotX;
  ngl::Mat4 rotY;

  // Set the initial model matrix and rotation matrix
  modelMatrix.identity();
  rotX.rotateX(m_spinXFace);
  rotY.rotateY(m_spinYFace);

  mouseGlobalTX = rotY*rotX;

  // Pass the projection, model, rotation and view matrices and light positions to the shader
  // First used by the bounding box
  shader->setRegisteredUniform("u_Projection", m_cam.getProjectionMatrix());
  shader->setRegisteredUniform("u_MV", modelMatrix * mouseGlobalTX * m_cam.getViewMatrix());

  shader->setRegisteredUniform("u_Light.Position", mouseGlobalTX * (m_cam.getEye() + ngl::Vec3(0.0f, 2.0f, 0.f)));
  shader->setRegisteredUniform("u_BackLight.Position", mouseGlobalTX * (m_cam.getEye() * ngl::Vec3(1.f, 1.f, -1.f) + ngl::Vec3(0.0f, 2.0f, 0.f)));

  // Execute the fluid system and simulation if simulation is enabled
  m_pbf.execute();

  // Loop through the particles and modify the model matrix to translate and scale the particles
  // to their respective locations and scales. Using a ngl::VAOPrimitive sphere to draw the particles
  std::vector<Particle *> particles = m_pbf.getParticles();
  for(unsigned int i = 0; i < particles.size(); ++i)
  {
    modelMatrix.identity();
    modelMatrix.scale(particles[i]->m_radius, particles[i]->m_radius, particles[i]->m_radius);
    modelMatrix.translate(particles[i]->m_pos.m_x, particles[i]->m_pos.m_y, particles[i]->m_pos.m_z);
    shader->setRegisteredUniform4f("u_Color", particles[i]->m_colour.m_x,  particles[i]->m_colour.m_y,  particles[i]->m_colour.m_z,  particles[i]->m_colour.m_w);
    shader->setRegisteredUniform("u_MV", modelMatrix * mouseGlobalTX * m_cam.getViewMatrix());
    particle->draw("particle");
  }

  // Record the frame end time
  m_end = std::chrono::system_clock::now();
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mouseMoveEvent (QMouseEvent * _event)
{
  // note the method buttons() is the button state when event was called
  // this is different from button() which is used to check which button was
  // pressed when the mousePress/Release event is generated
  if(m_rotate && _event->buttons() == Qt::LeftButton)
  {
    int diffx = _event->x() - m_origX;
    int diffy = _event->y() - m_origY;
    m_spinXFace += (float) 0.5f * diffy;
    m_spinYFace += (float) 0.5f * diffx;
    m_origX = _event->x();
    m_origY = _event->y();
  }
}


//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mousePressEvent ( QMouseEvent * _event)
{
  // this method is called when the mouse button is pressed in this case we
  // store the value where the maouse was clicked (x,y) and set the Rotate flag to true
  if(_event->button() == Qt::LeftButton)
  {
    m_origX = _event->x();
    m_origY = _event->y();
    m_rotate = true;
  }
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::mouseReleaseEvent ( QMouseEvent * _event )
{
  // this e  rotX++;vent is called when the mouse button is released
  // we then set Rotate to false
  if(_event->button() == Qt::LeftButton)
  {
    m_rotate = false;
  }
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::wheelEvent(QWheelEvent *_event)
{

}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::keyPressEvent(QKeyEvent *_event)
{
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
  {
    // escape key to quite
    case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
    // 1 to toggle simulation on/off
    case Qt::Key_1 : m_pbf.toggleSimulation(); break;
    // 2 to toggle wave machine on/off
    case Qt::Key_2 : m_pbf.toggleWaves(); break;
    default : break;
  }
  // finally update the GLWindow and re-draw

  update();
}
