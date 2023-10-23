/* DO NOT EDIT.
 * Generated by ../bin/vtkEncodeString
 * 
 * Define the vtkGPUVolumeRayCastMapper_MIPFourDependentNoCroppingFS string.
 *
 * Generated from file: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VolumeRendering/vtkGPUVolumeRayCastMapper_MIPFourDependentNoCroppingFS.glsl
 */
#include "vtkGPUVolumeRayCastMapper_MIPFourDependentNoCroppingFS.h"
const char *vtkGPUVolumeRayCastMapper_MIPFourDependentNoCroppingFS =
"/*=========================================================================\n"
"\n"
"  Program:   Visualization Toolkit\n"
"  Module:    vtkGPUVolumeRayCastMapper_MIPFourDependentNoCroppingFS.glsl\n"
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
"// Implementation of some functions used by the 4-component Maximum Intensity\n"
"// Projection (MIP) method when cropping is off.\n"
"\n"
"#version 110\n"
"\n"
"float initialMaxValue()\n"
"{\n"
"  return 0.0;\n"
"}\n"
"\n"
"vec4 initialColor()\n"
"{\n"
"  return vec4(0.0,0.0,0.0,0.0);\n"
"}\n"
"\n"
"void writeColorAndMaxScalar(vec4 color,\n"
"                            vec4 opacity,\n"
"                            float maxValue)\n"
"{\n"
"  // maxValue is not used\n"
"  \n"
"  // color framebuffer\n"
"  gl_FragColor.r = color.r*opacity.a;\n"
"  gl_FragColor.g = color.g*opacity.a;\n"
"  gl_FragColor.b = color.b*opacity.a;\n"
"  gl_FragColor.a=opacity.a;\n"
"}\n"
"\n";
