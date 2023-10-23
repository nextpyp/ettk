/* DO NOT EDIT.
 * Generated by ../bin/vtkEncodeString
 * 
 * Define the vtkGPUVolumeRayCastMapper_AdditiveNoCroppingFS string.
 *
 * Generated from file: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VolumeRendering/vtkGPUVolumeRayCastMapper_AdditiveNoCroppingFS.glsl
 */
#include "vtkGPUVolumeRayCastMapper_AdditiveNoCroppingFS.h"
const char *vtkGPUVolumeRayCastMapper_AdditiveNoCroppingFS =
"/*=========================================================================\n"
"\n"
"  Program:   Visualization Toolkit\n"
"  Module:    vtkGPUVolumeRayCastMapper_AdditiveNoCroppingFS.glsl\n"
"\n"
"  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen\n"
"  All rights reserved.\n"
"  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.\n"
"\n"
"     This software is distributed WITHOUT ANY WARRANTY; without even\n"
"     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR\n"
"     PURPOSE.  See the above copyright notice for more information.\n"
"\n"
"=========================================================================*/\n"
"\n"
"// Implementation of some functions used by the Additive method when cropping\n"
"// is off.\n"
"\n"
"#version 110\n"
"\n"
"float initialValue()\n"
"{\n"
"  return 0.0;\n"
"}\n"
"void writeColorAndSumScalar(vec4 color,\n"
"                            float sumValue)\n"
"{\n"
"  // we don't need to write sumValue to a buffer when there is no cropping.\n"
"  // color framebuffer\n"
"  gl_FragColor=color;\n"
"}\n"
"\n";
