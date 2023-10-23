// Loadable shader code
//
// Generated by ../../bin/ProcessShader
//
#ifndef __vtkShaderGLSL_h
#define __vtkShaderGLSL_h


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/TestAppVarFrag.glsl
static const char* vtkShaderGLSLTestAppVarFragCode0 =
"//\n"
"// Begin \"3Dlabs-License.txt\"\n"
"//\n"
"// Copyright (C) 2002-2005  3Dlabs Inc. Ltd.\n"
"// All rights reserved.\n"
"//\n"
"// Redistribution and use in source and binary forms, with or without\n"
"// modification, are permitted provided that the following conditions\n"
"// are met:\n"
"//\n"
"//     Redistributions of source code must retain the above copyright\n"
"//     notice, this list of conditions and the following disclaimer.\n"
"//\n"
"//     Redistributions in binary form must reproduce the above\n"
"//     copyright notice, this list of conditions and the following\n"
"//     disclaimer in the documentation and/or other materials provided\n"
"//     with the distribution.\n"
"//\n"
"//     Neither the name of 3Dlabs Inc. Ltd. nor the names of its\n"
"//     contributors may be used to endorse or promote products derived\n"
"//     from this software without specific prior written permission.\n"
"//\n"
"// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"// \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
"// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
"// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
"// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
"// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
"// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN\n"
"// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"// POSSIBILITY OF SUCH DAMAGE.\n"
"//\n"
"//\n"
"// End \"3Dlabs-License.txt\"\n"
"\n"
"\n"
"\n"
"uniform vec3 Color;\n"
"uniform vec3 AmbientColor;\n"
"uniform vec3 DiffuseColor;\n"
"uniform vec3 SpecularColor;\n"
"uniform vec3 EdgeColor;\n"
"\n"
"uniform float Ambient;\n"
"uniform float Diffuse;\n"
"uniform float Specular;\n"
"uniform float SpecularPower;\n"
"uniform float Opacity;\n"
"\n"
"uniform float PointSize;\n"
"uniform float LineWidth;\n"
"\n"
"uniform int LineStipplePattern;\n"
"uniform int LineStippleRepeatFactor;\n"
"uniform int Interpolation;\n"
"uniform int Representation;\n"
"uniform int EdgeVisibility;\n"
"uniform int BackfaceCulling;\n"
"uniform int FrontfaceCulling;\n"
"\n"
"\n"
"uniform vec3  SurfaceColor; // (0.75, 0.75, 0.75)\n"
"uniform vec3  WarmColor;    // (0.6, 0.6, 0.0)\n"
"uniform vec3  CoolColor;    // (0.0, 0.0, 0.6)\n"
"uniform float DiffuseWarm;  // 0.45\n"
"uniform float DiffuseCool;  // 0.45\n"
"\n"
"varying float NdotL;\n"
"varying vec3  ReflectVec;\n"
"varying vec3  ViewVec;\n"
"\n"
"uniform vec4 appVara;\n"
"uniform vec4 appVarb;\n"
"uniform vec4 appVarc;\n"
"uniform vec4 appVard;\n"
"uniform vec4 appVare;\n"
"uniform vec4 appVarf;\n"
"uniform vec4 appVarg;\n"
"\n"
"void main (void)\n"
"{\n"
"    vec3 kcool    = min(CoolColor + DiffuseCool * SurfaceColor, 1.0);\n"
"    vec3 kwarm    = min(WarmColor + DiffuseWarm * SurfaceColor, 1.0); \n"
"    vec3 kfinal   = mix(kcool, kwarm, NdotL);\n"
"\n"
"    vec3 nreflect = normalize(ReflectVec);\n"
"    vec3 nview    = normalize(ViewVec);\n"
"\n"
"    float spec    = max(dot(nreflect, nview), 0.0);\n"
"    spec          = pow(spec, 32.0);\n"
"\n"
"    gl_FragColor = vec4 (min(kfinal + spec, 1.0), 1.0);\n"
"\n"
"\n"
"    if( 0\n"
"        || appVara.x != 0.37714 || appVara.y != 0.61465 || appVara.z != 0.48399 || appVara.w != 0.68252\n"
"        || appVarb.x != 0.03900 || appVarb.y != 0.15857 || appVarb.z != 0.57913 || appVarb.w != 0.54458\n"
"        || appVarc.x != 0.97061 || appVarc.y != 0.86053 || appVarc.z != 0.63583 || appVarc.w != 0.51058\n"
"        || appVard.x != 0.12885 || appVard.y != 0.91490 || appVard.z != 0.86394 || appVard.w != 0.58951\n"
"        || appVare.x != 0.23403 || appVare.y != 0.35340 || appVare.z != 0.52559 || appVare.w != 0.77830\n"
"        || appVarf.x != 0.19550 || appVarf.y != 0.17429 || appVarf.z != 0.89958 || appVarf.w != 0.15063\n"
"        || appVarg.x != 0.75796 || appVarg.y != 0.48072 || appVarg.z != 0.07728 || appVarg.w != 0.16434\n"
"      )\n"
"      {\n"
"      gl_FragColor = vec4 (1.0, 0.0, 0.0, 1.0);\n"
"      }\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTestAppVarFragGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTestAppVarFragCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTestAppVarFragCode0);
  return res;
}


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/TestVertex.glsl
static const char* vtkShaderGLSLTestVertexCode0 =
"//\n"
"// Begin \"3Dlabs-License.txt\"\n"
"//\n"
"// Copyright (C) 2002-2005  3Dlabs Inc. Ltd.\n"
"// All rights reserved.\n"
"//\n"
"// Redistribution and use in source and binary forms, with or without\n"
"// modification, are permitted provided that the following conditions\n"
"// are met:\n"
"//\n"
"//     Redistributions of source code must retain the above copyright\n"
"//     notice, this list of conditions and the following disclaimer.\n"
"//\n"
"//     Redistributions in binary form must reproduce the above\n"
"//     copyright notice, this list of conditions and the following\n"
"//     disclaimer in the documentation and/or other materials provided\n"
"//     with the distribution.\n"
"//\n"
"//     Neither the name of 3Dlabs Inc. Ltd. nor the names of its\n"
"//     contributors may be used to endorse or promote products derived\n"
"//     from this software without specific prior written permission.\n"
"//\n"
"// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"// \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
"// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
"// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
"// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
"// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
"// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN\n"
"// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"// POSSIBILITY OF SUCH DAMAGE.\n"
"//\n"
"//\n"
"// End \"3Dlabs-License.txt\"\n"
"\n"
"uniform vec3  LightPosition;  // (0.0, 10.0, 4.0) \n"
"\n"
"varying float NdotL;\n"
"varying vec3  ReflectVec;\n"
"varying vec3  ViewVec;\n"
"\n"
"struct light\n"
"{\n"
"  uniform vec3 position;\n"
"  uniform vec3 color;\n"
"};\n"
"\n"
"void main(void)\n"
"{\n"
"  uniform light l1;\n"
"  uniform light l2;\n"
"    vec3 ecPos      = vec3 (gl_ModelViewMatrix * gl_Vertex);\n"
"    vec3 tnorm      = normalize(gl_NormalMatrix * gl_Normal);\n"
"    vec3 lightVec   = normalize(LightPosition - ecPos);\n"
"    ReflectVec      = normalize(reflect(-lightVec, tnorm));\n"
"    ViewVec         = normalize(-ecPos);\n"
"    NdotL           = (dot(lightVec, tnorm) + 1.0) * 0.5;\n"
"    gl_Position     = ftransform();\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTestVertexGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTestVertexCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTestVertexCode0);
  return res;
}


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/TestVtkPropertyFrag.glsl
static const char* vtkShaderGLSLTestVtkPropertyFragCode0 =
"//\n"
"// Begin \"3Dlabs-License.txt\"\n"
"//\n"
"// Copyright (C) 2002-2005  3Dlabs Inc. Ltd.\n"
"// All rights reserved.\n"
"//\n"
"// Redistribution and use in source and binary forms, with or without\n"
"// modification, are permitted provided that the following conditions\n"
"// are met:\n"
"//\n"
"//     Redistributions of source code must retain the above copyright\n"
"//     notice, this list of conditions and the following disclaimer.\n"
"//\n"
"//     Redistributions in binary form must reproduce the above\n"
"//     copyright notice, this list of conditions and the following\n"
"//     disclaimer in the documentation and/or other materials provided\n"
"//     with the distribution.\n"
"//\n"
"//     Neither the name of 3Dlabs Inc. Ltd. nor the names of its\n"
"//     contributors may be used to endorse or promote products derived\n"
"//     from this software without specific prior written permission.\n"
"//\n"
"// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"// \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
"// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
"// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
"// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
"// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
"// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN\n"
"// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"// POSSIBILITY OF SUCH DAMAGE.\n"
"//\n"
"//\n"
"// End \"3Dlabs-License.txt\"\n"
"\n"
"\n"
"\n"
"uniform vec3 Color;\n"
"uniform vec3 AmbientColor;\n"
"uniform vec3 DiffuseColor;\n"
"uniform vec3 SpecularColor;\n"
"uniform vec3 EdgeColor;\n"
"\n"
"uniform float Ambient;\n"
"uniform float Diffuse;\n"
"uniform float Specular;\n"
"uniform float SpecularPower;\n"
"uniform float Opacity;\n"
"\n"
"uniform float PointSize;\n"
"uniform float LineWidth;\n"
"\n"
"uniform int LineStipplePattern;\n"
"uniform int LineStippleRepeatFactor;\n"
"uniform int Interpolation;\n"
"uniform int Representation;\n"
"uniform int EdgeVisibility;\n"
"uniform int BackfaceCulling;\n"
"uniform int FrontfaceCulling;\n"
"\n"
"\n"
"uniform vec3  SurfaceColor; // (0.75, 0.75, 0.75)\n"
"uniform vec3  WarmColor;    // (0.6, 0.6, 0.0)\n"
"uniform vec3  CoolColor;    // (0.0, 0.0, 0.6)\n"
"uniform float DiffuseWarm;  // 0.45\n"
"uniform float DiffuseCool;  // 0.45\n"
"\n"
"varying float NdotL;\n"
"varying vec3  ReflectVec;\n"
"varying vec3  ViewVec;\n"
"\n"
"void main (void)\n"
"{\n"
"    vec3 kcool    = min(CoolColor + DiffuseCool * SurfaceColor, 1.0);\n"
"    vec3 kwarm    = min(WarmColor + DiffuseWarm * SurfaceColor, 1.0); \n"
"    vec3 kfinal   = mix(kcool, kwarm, NdotL);\n"
"\n"
"    vec3 nreflect = normalize(ReflectVec);\n"
"    vec3 nview    = normalize(ViewVec);\n"
"\n"
"    float spec    = max(dot(nreflect, nview), 0.0);\n"
"    spec          = pow(spec, 32.0);\n"
"\n"
"    gl_FragColor = vec4 (min(kfinal + spec, 1.0), 1.0);\n"
"\n"
"\n"
"    if( 0\n"
"      || AmbientColor.x!=0.75 || AmbientColor.y!=0.751 || AmbientColor.z!=0.752\n"
"      || DiffuseColor.x!=0.61 || DiffuseColor.y!=0.62 || DiffuseColor.z!=0.006\n"
"      || SpecularColor.x!=0.001 || SpecularColor.y!=0.002 || SpecularColor.z!=0.61\n"
"      || EdgeColor.x!=0.1 || EdgeColor.y!=0.2 || EdgeColor.z!=0.3\n"
"      || Ambient!=0.45\n"
"      || Diffuse!=0.451\n"
"      || Specular!=0.4\n"
"      || SpecularPower!=1.0\n"
"      || Opacity!=1.0\n"
"      || PointSize!=1.0\n"
"      || LineWidth!=1.0\n"
"      || LineStipplePattern!=0\n"
"      || LineStippleRepeatFactor!=1\n"
"      || Interpolation!=1\n"
"      || Representation!=2\n"
"      || EdgeVisibility!=0\n"
"      || BackfaceCulling!=0\n"
"      || FrontfaceCulling!=0\n"
"    )\n"
"      {\n"
"      gl_FragColor = vec4 (1.0, 0.0, 0.0, 1.0);\n"
"      }\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTestVtkPropertyFragGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTestVtkPropertyFragCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTestVtkPropertyFragCode0);
  return res;
}


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/TestMatrixFrag.glsl
static const char* vtkShaderGLSLTestMatrixFragCode0 =
"//\n"
"// Begin \"3Dlabs-License.txt\"\n"
"//\n"
"// Copyright (C) 2002-2005  3Dlabs Inc. Ltd.\n"
"// All rights reserved.\n"
"//\n"
"// Redistribution and use in source and binary forms, with or without\n"
"// modification, are permitted provided that the following conditions\n"
"// are met:\n"
"//\n"
"//     Redistributions of source code must retain the above copyright\n"
"//     notice, this list of conditions and the following disclaimer.\n"
"//\n"
"//     Redistributions in binary form must reproduce the above\n"
"//     copyright notice, this list of conditions and the following\n"
"//     disclaimer in the documentation and/or other materials provided\n"
"//     with the distribution.\n"
"//\n"
"//     Neither the name of 3Dlabs Inc. Ltd. nor the names of its\n"
"//     contributors may be used to endorse or promote products derived\n"
"//     from this software without specific prior written permission.\n"
"//\n"
"// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"// \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
"// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
"// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
"// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
"// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
"// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN\n"
"// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"// POSSIBILITY OF SUCH DAMAGE.\n"
"//\n"
"//\n"
"// End \"3Dlabs-License.txt\"\n"
"\n"
"uniform mat2 testMat2;\n"
"uniform mat3 testMat3;\n"
"uniform mat4 testMat4;\n"
"\n"
"\n"
"\n"
"\n"
"uniform vec3  SurfaceColor; // (0.75, 0.75, 0.75)\n"
"uniform vec3  WarmColor;    // (0.6, 0.6, 0.0)\n"
"uniform vec3  CoolColor;    // (0.0, 0.0, 0.6)\n"
"uniform float DiffuseWarm;  // 0.45\n"
"uniform float DiffuseCool;  // 0.45\n"
"\n"
"varying float NdotL;\n"
"varying vec3  ReflectVec;\n"
"varying vec3  ViewVec;\n"
"\n"
"void main (void)\n"
"{\n"
"    vec3 kcool    = min(CoolColor + DiffuseCool * SurfaceColor, 1.0);\n"
"    vec3 kwarm    = min(WarmColor + DiffuseWarm * SurfaceColor, 1.0); \n"
"    vec3 kfinal   = mix(kcool, kwarm, NdotL);\n"
"\n"
"    vec3 nreflect = normalize(ReflectVec);\n"
"    vec3 nview    = normalize(ViewVec);\n"
"\n"
"    float spec    = max(dot(nreflect, nview), 0.0);\n"
"    spec          = pow(spec, 32.0);\n"
"\n"
"    gl_FragColor = vec4 (min(kfinal + spec, 1.0), 1.0);\n"
"\n"
"\n"
"    if( 0\n"
"\n"
"      || testMat2[0][0]!=3.14159 || testMat2[0][1]!=6.28319\n"
"      || testMat2[1][0]!=9.42478 || testMat2[1][1]!=12.5664\n"
"\n"
"      || testMat3[0][0]!=9.4248  || testMat3[0][1]!=18.8496 || testMat3[0][2]!=28.2743\n"
"      || testMat3[1][0]!=37.6991 || testMat3[1][1]!=47.1239 || testMat3[1][2]!=56.5487\n"
"      || testMat3[2][0]!=65.9734 || testMat3[2][1]!=75.3982 || testMat3[2][2]!=84.8230\n"
"\n"
"      || testMat4[0][0]!=12.5664  || testMat4[0][1]!=25.1327  || testMat4[0][2]!=37.6991  || testMat4[0][3]!=50.2655 \n"
"      || testMat4[1][0]!=62.8319  || testMat4[1][1]!=75.3982  || testMat4[1][2]!=87.9646  || testMat4[1][3]!=100.5310\n"
"      || testMat4[2][0]!=113.0973 || testMat4[2][1]!=125.6637 || testMat4[2][2]!=138.2301 || testMat4[2][3]!=150.7964\n"
"      || testMat4[3][0]!=163.3628 || testMat4[3][1]!=175.9292 || testMat4[3][2]!=188.4956 || testMat4[3][3]!=201.0619\n"
"\n"
"    )\n"
"      {\n"
"      gl_FragColor = vec4 (1.0, 0.0, 0.0, 1.0);\n"
"      }\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTestMatrixFragGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTestMatrixFragCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTestMatrixFragCode0);
  return res;
}


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/TestScalarVectorFrag.glsl
static const char* vtkShaderGLSLTestScalarVectorFragCode0 =
"//\n"
"// Begin \"3Dlabs-License.txt\"\n"
"//\n"
"// Copyright (C) 2002-2005  3Dlabs Inc. Ltd.\n"
"// All rights reserved.\n"
"//\n"
"// Redistribution and use in source and binary forms, with or without\n"
"// modification, are permitted provided that the following conditions\n"
"// are met:\n"
"//\n"
"//     Redistributions of source code must retain the above copyright\n"
"//     notice, this list of conditions and the following disclaimer.\n"
"//\n"
"//     Redistributions in binary form must reproduce the above\n"
"//     copyright notice, this list of conditions and the following\n"
"//     disclaimer in the documentation and/or other materials provided\n"
"//     with the distribution.\n"
"//\n"
"//     Neither the name of 3Dlabs Inc. Ltd. nor the names of its\n"
"//     contributors may be used to endorse or promote products derived\n"
"//     from this software without specific prior written permission.\n"
"//\n"
"// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"// \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
"// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
"// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
"// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
"// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
"// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN\n"
"// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"// POSSIBILITY OF SUCH DAMAGE.\n"
"//\n"
"//\n"
"// End \"3Dlabs-License.txt\"\n"
"\n"
"uniform float testFloat;\n"
"uniform vec2  testVec2;\n"
"uniform vec3  testVec3;\n"
"uniform vec4  testVec4;\n"
"\n"
"uniform int   testInt;\n"
"uniform ivec2  testIVec2;\n"
"uniform ivec3  testIVec3;\n"
"uniform ivec4  testIVec4;\n"
"\n"
"uniform mat2 testMat2;\n"
"uniform mat3 testMat3;\n"
"uniform mat4 testMat4;\n"
"\n"
"struct tStruct {\n"
"  float f;\n"
"  vec2 f2;\n"
"  vec3 f3;\n"
"  vec4 f4;\n"
"};\n"
"\n"
"uniform tStruct tStruct2;\n"
"\n"
"\n"
"uniform vec3  SurfaceColor; // (0.75, 0.75, 0.75)\n"
"uniform vec3  WarmColor;    // (0.6, 0.6, 0.0)\n"
"uniform vec3  CoolColor;    // (0.0, 0.0, 0.6)\n"
"uniform float DiffuseWarm;  // 0.45\n"
"uniform float DiffuseCool;  // 0.45\n"
"\n"
"varying float NdotL;\n"
"varying vec3  ReflectVec;\n"
"varying vec3  ViewVec;\n"
"\n"
"void main (void)\n"
"{\n"
"    vec3 kcool    = min(CoolColor + DiffuseCool * SurfaceColor, 1.0);\n"
"    vec3 kwarm    = min(WarmColor + DiffuseWarm * SurfaceColor, 1.0); \n"
"    vec3 kfinal   = mix(kcool, kwarm, NdotL);\n"
"\n"
"    vec3 nreflect = normalize(ReflectVec);\n"
"    vec3 nview    = normalize(ViewVec);\n"
"\n"
"    float spec    = max(dot(nreflect, nview), 0.0);\n"
"    spec          = pow(spec, 32.0);\n"
"\n"
"    gl_FragColor = vec4 (min(kfinal + spec, 1.0), 1.0);\n"
"\n"
"\n"
"    if( 0\n"
"      || testFloat!=1.0\n"
"      || testVec2.x!=1.0 || testVec2.y!=2.0\n"
"      || testVec3.x!=1.0 || testVec3.y!=2.0 || testVec3.z!=3.0\n"
"      || testVec4.x!=1.0 || testVec4.y!=2.0 || testVec4.z!=3.0 || testVec4.w!=4.0\n"
"\n"
"      || testInt!=1\n"
"      || testIVec2.x!=1 || testIVec2.y!=2\n"
"      || testIVec3.x!=1 || testIVec3.y!=2 || testIVec3.z!=3\n"
"      || testIVec4.x!=1 || testIVec4.y!=2 || testIVec4.z!=3 || testIVec4.w!=4\n"
"\n"
"      || testMat2[0][0]!=3.14159 || testMat2[0][1]!=6.28319\n"
"      || testMat2[1][0]!=9.42478 || testMat2[1][1]!=12.5664\n"
"\n"
"      || testMat3[0][0]!=9.4248  || testMat3[0][1]!=18.8496 || testMat3[0][2]!=28.2743\n"
"      || testMat3[1][0]!=37.6991 || testMat3[1][1]!=47.1239 || testMat3[1][2]!=56.5487\n"
"      || testMat3[2][0]!=65.9734 || testMat3[2][1]!=75.3982 || testMat3[2][2]!=84.8230\n"
"\n"
"      || testMat4[0][0]!=12.5664  || testMat4[0][1]!=25.1327  || testMat4[0][2]!=37.6991  || testMat4[0][3]!=50.2655\n"
"      || testMat4[1][0]!=62.8319  || testMat4[1][1]!=75.3982  || testMat4[1][2]!=87.9646  || testMat4[1][3]!=100.5310\n"
"      || testMat4[2][0]!=113.0973 || testMat4[2][1]!=125.6637 || testMat4[2][2]!=138.2301 || testMat4[2][3]!=150.7964\n"
"      || testMat4[3][0]!=163.3628 || testMat4[3][1]!=175.9292 || testMat4[3][2]!=188.4956 || testMat4[3][3]!=201.0619\n"
"\n"
"      || tStruct2.f!=2.0\n"
"      || tStruct2.f2.x!=2.0 || tStruct2.f2.y!=2.1\n"
"      || tStruct2.f3.x!=2.0 || tStruct2.f3.y!=2.1 || tStruct2.f3.z!=2.2\n"
"      || tStruct2.f4.x!=2.0 || tStruct2.f4.y!=2.1 || tStruct2.f4.z!=2.2 || tStruct2.f4.w!=2.3\n"
"\n"
"    )\n"
"      {\n"
"      gl_FragColor = vec4 (1.0, 0.0, 0.0, 1.0);\n"
"      }\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTestScalarVectorFragGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTestScalarVectorFragCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTestScalarVectorFragCode0);
  return res;
}


// From file /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary/GLSLShaders/Twisted.glsl
static const char* vtkShaderGLSLTwistedCode0 =
"uniform float Rate;\n"
"\n"
"void main()\n"
"{\n"
"  float angle = gl_Vertex[2];\n"
"  float rate = Rate;\n"
"  vec4 newPos;\n"
"  newPos[0] = cos(angle* rate) * gl_Vertex[0] + sin(angle* rate) *gl_Vertex[1];\n"
"  newPos[1] = -sin(angle* rate) * gl_Vertex[0] + cos(angle* rate) *gl_Vertex[1];\n"
"  newPos[2] = gl_Vertex[2];\n"
"  newPos[3] = gl_Vertex[3];\n"
"\n"
"  vec3 newNormal;\n"
"  newNormal[0] = cos(angle* rate) * gl_Normal[0] + sin(angle* rate) *gl_Normal[1];\n"
"  newNormal[1] = -sin(angle* rate) * gl_Normal[0] + cos(angle* rate) *gl_Normal[1];\n"
"  newNormal[2] = gl_Normal[2];\n"
"\n"
"  vec4 ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;\n"
"  vec3 normal = normalize(gl_NormalMatrix * newNormal);\n"
"  vec3 lightDir = normalize(vec3(gl_LightSource[0].position));\n"
"  float NdotL = max(dot(normal, lightDir), 0.0);\n"
"\n"
"  vec4 col = vec4(0,0.5,0.5,1);\n"
"  vec4 diffuse = col * gl_LightSource[0].diffuse * NdotL;\n"
"  vec4 specular;\n"
"  \n"
"  if (NdotL > 0.0) \n"
"    {\n"
"    float NdotHV = max(dot(normal, gl_LightSource[0].halfVector.xyz), 0.0);\n"
"    specular = gl_FrontMaterial.specular * gl_LightSource[0].specular * \n"
"      pow(NdotHV, gl_FrontMaterial.shininess);\n"
"    }\n"
"  gl_FrontColor = diffuse + ambient + specular;\n"
"\n"
"  gl_Position = gl_ModelViewProjectionMatrix * newPos;\n"
"\n"
"\n"
"\n"
"}\n"
"\n";
// Get single string
char* vtkShaderGLSLTwistedGetCode()
{
  int len = ( 0
    + static_cast<int>(strlen(vtkShaderGLSLTwistedCode0)) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, vtkShaderGLSLTwistedCode0);
  return res;
}



#endif
