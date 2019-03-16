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

#ifndef _COMMON_HEADERS_H_
#define _COMMON_HEADERS_H_

// precision
#define HIGH_PRECISION

// matlab debugger
// NOTE: if you don't have Matlab R2010a or higher versions installed on your computer, please comment the next line.
// #define ENABLE_MATLAB_DEBUGGING

// single or double presicion
#ifdef HIGH_PRECISION
    typedef double ScalarType;
    #define TW_TYPE_SCALAR_TYPE TW_TYPE_DOUBLE
#else
    typedef float ScalarType;
    #define TW_TYPE_SCALAR_TYPE TW_TYPE_FLOAT
#endif

// small number and large number
#ifdef HIGH_PRECISION
    #define EPSILON 1e-15
#else
    #define EPSILON 1e-6
#endif

#define LARGER_EPSILON 1e-6

// default values
#define DEFAULT_SCREEN_WIDTH 1024
#define DEFAULT_SCREEN_HEIGHT 768

// selection radius
#define DEFAULT_SELECTION_RADIUS 0.1

// file localtions
#define DEFAULT_VERT_SHADER_FILE "./shaders/vert.glsl"
#define DEFAULT_FRAG_SHADER_FILE "./shaders/frag.glsl"
#define DEFAULT_SCENE_FILE "./scenes/test_scene.xml"
#define DEFAULT_CONFIG_FILE "./config/config.txt"
#define DEFAULT_CLOTH_STATE_FILE "./config/cloth.txt"
#define DEFAULT_TEXTURE_FILE "./textures/chinese_porcelain_1.png"
#define DEFAULT_SCREEN_SHOT_FILE "./output/ScreenShot.png"
#define DEFAULT_OUTPUT_CLOTH_STATE_FILE "./output/cloth.txt"
#define DEFAULT_OUTPUT_OBJ_FILE "./output/cloth.obj"
#define DEFAULT_MODEL "./mesh_models/frog_reg.mesh"
#define DEFAULT_OBJ_MODEL "./obj_models/bunny.obj"

#endif