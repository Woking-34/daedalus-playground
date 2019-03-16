// ---------------------------------------------------------------------------------//
// Copyright (c) 2013, Regents of the University of Pennsylvania                    //
// All rights reserved.                                                             //
//                                                                                  //
// Redistribution and use in source and binary forms, with or without               //
// modification, are permitted provided that the following conditions are met:      //
//     * Redistributions of source code must retain the above copyright             //
//       notice, this list of conditions and the following disclaimer.              //
//     * Redistributions in binary form must reproduce the above copyright          //
//       notice, this list of conditions and the following disclaimer in the        //
//       documentation and/or other materials provided with the distribution.       //
//     * Neither the name of the <organization> nor the                             //
//       names of its contributors may be used to endorse or promote products       //
//       derived from this software without specific prior written permission.      //
//                                                                                  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  //
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           //
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY               //
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       //
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     //
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      //
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       //
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//                                                                                  //
// Contact Tiantian Liu (ltt1598@gmail.com) if you have any questions.              //
//----------------------------------------------------------------------------------//

#include "camera.h"

Camera::Camera(void)
{
}

Camera::~Camera(void)
{

}

void Camera::Reset(int width, int height)
{
    // setup default camera parameters
    m_eye_distance = 15.0;
    m_head = 30.0;
    m_pitch = 45.0;
    m_lookat = glm::vec3(0.0, 4.5, 0.0);
    m_up = glm::vec3(0.0, 1.0, 0.0);
    m_fovy = 60.0;
    m_width = width;
    m_height = height;
    m_znear = 0.1;
    m_zfar = 500.0;

    updateViewMatrix();
    updateProjectionMatrix();
}

void Camera::Lookat(Mesh* mesh)
{
    unsigned int mid_index = mesh->m_vertices_number/2;

    EigenVector3 lookat = mesh->m_current_positions.block_vector(mid_index);
    m_lookat = glm::vec3(lookat[0], lookat[1], lookat[2]);
    updateViewMatrix();
}

void Camera::DrawAxis()
{
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH);

    // store previous states
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    // change view matrix.
    glm::vec3 axis_cam_pos = float(2.0) * glm::normalize(m_position - m_lookat);
    glLoadMatrixf(&(glm::lookAt(axis_cam_pos, glm::vec3(0.0, 0.0, 0.0), m_up)[0][0]));

    // change viewport
    glViewport(m_width * 15 / 16, 0, m_width / 16, m_height / 16);

    //Draw axis.
    glBegin(GL_LINES);
    glColor3d(1.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(1.0, 0.0, 0.0);

    glColor3d(0.0, 1.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 1.0, 0.0);

    glColor3d(0.0, 0.0, 1.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 1.0);
    glEnd();

    // restore everything
    glViewport(0, 0, m_width, m_height);
    glPopMatrix();
    glPopAttrib();
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH);
}

// mouse interactions
void Camera::MouseChangeDistance(float coe, float dx, float dy)
{
    m_eye_distance -= dy * coe;
    if (m_eye_distance < 4.0) m_eye_distance = 4.0; 
    updateViewMatrix();
}
void Camera::MouseChangeLookat(float coe, float dx, float dy)
{
    glm::vec3 vdir(m_lookat - m_position);
    glm::vec3 u(glm::normalize(glm::cross(vdir, m_up)));
    glm::vec3 v(glm::normalize(glm::cross(u, vdir)));

    m_lookat += coe * (dy * v - dx * u);
    updateViewMatrix();
}
void Camera::MouseChangeHeadPitch(float coe, float dx, float dy)
{
    m_head += dy * coe;
    m_pitch += dx * coe;

    updateViewMatrix();
}

// resize
void Camera::ResizeWindow(int w, int h)
{
    this->m_width = w;
    this->m_height = h;

    updateProjectionMatrix();
}

glm::vec3 Camera::GetRaycastDirection(int mouse_x, int mouse_y)
{
    float x = (float)(mouse_x) / (float)(m_width-1); 
    float y = 1.0 - (float)(mouse_y) / (float)(m_height-1); 

    // Viewing vector
    glm::vec3 E = m_position;
    glm::vec3 U = m_up;
    glm::vec3 C = glm::normalize(m_lookat - m_position); // implies viewing plane distancei s 1

    float phi = glm::radians(m_fovy/2.0);

    // Vector A = C x U
    glm::vec3 A = glm::normalize(glm::cross(C, U));
    // The REAL up vector B = A x C
    glm::vec3 B = glm::normalize(glm::cross(A, C));
    // View Center M = E + C
    glm::vec3 M = E + C;

    // V || B, but on NCD
    glm::vec3 V = B * glm::tan(phi);
    // H || A, but on NDC
    // If you didn't use theta here, you can simply use the ratio between this->width() and this->height()
    glm::vec3 H = A * glm::tan(phi) / (float)m_height * (float)m_width;

    // Clicking point on the screen. World Coordinate.
    glm::vec3 P = M + float(2.0*x - 1.0)*H + float(2.0*y - 1.0)*V;

    m_cached_projection_plane_center = M;
    m_cached_projection_plane_xdir = H;
    m_cached_projection_plane_ydir = V;

    glm::vec3 dir = glm::normalize(P-E);

    return dir;
}

glm::vec3 Camera::GetCurrentTargetPoint(int mouse_x, int mouse_y)
{
    // assume camera is not moving
    float x = (float)(mouse_x) / (float)(m_width-1); 
    float y = 1.0f - (float)(mouse_y) / (float)(m_height-1); 

    glm::vec3 P = m_cached_projection_plane_center + float(2.0f*x - 1.0f)*m_cached_projection_plane_xdir + float(2.0f*y - 1.0f)*m_cached_projection_plane_ydir;

    glm::vec3 dir = P-m_position;

    glm::vec3 fixed_point_glm = m_cached_projection_plane_distance*dir+m_position;
    
    return fixed_point_glm;
}

// private field
void Camera::updateViewMatrix()
{
    float r_head = glm::radians(m_head), r_pitch = glm::radians(m_pitch);
    m_position.x = m_lookat.x + m_eye_distance * glm::cos(r_head) * glm::cos(r_pitch);
    m_position.y = m_lookat.y + m_eye_distance * glm::sin(r_head);
    m_position.z = m_lookat.z + m_eye_distance * glm::cos(r_head) * glm::sin(r_pitch);

    m_up = glm::vec3(0.0, (glm::cos(r_head) > 0.0) ? 1.0 : -1.0, 0.0);
    m_view = glm::lookAt(m_position, m_lookat, m_up);
}
void Camera::updateProjectionMatrix()
{
    m_projection = glm::perspective(m_fovy, static_cast<float>(m_width) / static_cast<float>(m_height), m_znear, m_zfar);
}