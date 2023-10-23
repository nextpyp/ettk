#include <ctype.h>  /* NOLINT */
#include <stdio.h>  /* NOLINT */
#include <stdlib.h> /* NOLINT */
#include <string.h> /* NOLINT */

#if defined(_MSC_VER)
#pragma warning(disable : 4996) /* deprecation */
#endif

#include "vtkTestDriver.h"


/* Forward declare test functions. */
int otherCoordinate(int, char*[]);
int TestPriorityStreaming(int, char*[]);
int LoadOpenGLExtension(int, char*[]);
int TestActorLightingFlag(int, char*[]);
int TestAnimationScene(int, char*[]);
int TestBackfaceCulling(int, char*[]);
int TestBlurAndSobelPasses(int, char*[]);
int TestDynamic2DLabelMapper(int, char*[]);
int TestFBO(int, char*[]);
int TestFollowerPicking(int, char*[]);
int TestGaussianBlurPass(int, char*[]);
int TestGlyph3DMapper(int, char*[]);
int TestGlyph3DMapperMasking(int, char*[]);
int TestGlyph3DMapperOrientationArray(int, char*[]);
int TestGlyph3DMapperPicking(int, char*[]);
int TestGPUInfo(int, char*[]);
int TestGradientBackground(int, char*[]);
int TestHomogeneousTransformOfActor(int, char*[]);
int TestImageResliceMapperAlpha(int, char*[]);
int TestImageResliceMapperBackground(int, char*[]);
int TestImageResliceMapperBorder(int, char*[]);
int TestImageResliceMapperInterpolation(int, char*[]);
int TestImageResliceMapperOffAxis(int, char*[]);
int TestImageResliceMapperOrient3D(int, char*[]);
int TestImageResliceMapperSlab(int, char*[]);
int TestImageSliceMapperAlpha(int, char*[]);
int TestImageSliceMapperBackground(int, char*[]);
int TestImageSliceMapperBorder(int, char*[]);
int TestImageSliceMapperOrient2D(int, char*[]);
int TestImageSliceMapperOrient3D(int, char*[]);
int TestImageSliceMapperInterpolation(int, char*[]);
int TestImageStack(int, char*[]);
int TestInteractorTimers(int, char*[]);
int TestLabelPlacer(int, char*[]);
int TestLabelPlacer2D(int, char*[]);
int TestLabelPlacerCoincidentPoints(int, char*[]);
int TestLabelPlacementMapper2D(int, char*[]);
int TestLabelPlacementMapperCoincidentPoints(int, char*[]);
int TestLightActor(int, char*[]);
int TestLODActor(int, char*[]);
int TestManyActors(int, char*[]);
int TestOrderedTriangulator(int, char*[]);
int TestOpacity(int, char*[]);
int TestOpenGLPolyDataMapper(int, char*[]);
int TestOSConeCxx(int, char*[]);
int TestPOVExporter(int, char*[]);
int TestResetCameraVerticalAspectRatio(int, char*[]);
int TestResetCameraVerticalAspectRatioParallel(int, char*[]);
int TestSetImageOrientation(int, char*[]);
int TestSobelGradientMagnitudePass(int, char*[]);
int TestShadowMapPass(int, char*[]);
int TestTextActorAlphaBlending(int, char*[]);
int TestTextActorDepthPeeling(int, char*[]);
int TestTextActor3DAlphaBlending(int, char*[]);
int TestTextActor3DDepthPeeling(int, char*[]);
int TestTexturedBackground(int, char*[]);
int TestTextureSize(int, char*[]);
int TestTDx(int, char*[]);
int TestTilingCxx(int, char*[]);
int TestTransformCoordinateUseDouble(int, char*[]);
int TestTranslucentLUTAlphaBlending(int, char*[]);
int TestTranslucentLUTDepthPeeling(int, char*[]);
int TestTranslucentLUTDepthPeelingPass(int, char*[]);
int TestTranslucentLUTTextureAlphaBlending(int, char*[]);
int TestTranslucentLUTTextureDepthPeeling(int, char*[]);
int TestGenericVertexAttributesGLSLCxx(int, char*[]);
int TestGenericVertexAttributesGLSLAlphaBlending(int, char*[]);
int TestGenericVertexAttributesGLSLDepthPeelingPass(int, char*[]);


#ifdef __cplusplus
#define CM_CAST(TYPE, EXPR) static_cast<TYPE>(EXPR)
#else
#define CM_CAST(TYPE, EXPR) (TYPE)(EXPR)
#endif

/* Create map.  */

typedef int (*MainFuncPointer)(int, char* []);
typedef struct
{
  const char* name;
  MainFuncPointer func;
} functionMapEntry;

static functionMapEntry cmakeGeneratedFunctionMapEntries[] = {
    {
    "otherCoordinate",
    otherCoordinate
  },
  {
    "TestPriorityStreaming",
    TestPriorityStreaming
  },
  {
    "LoadOpenGLExtension",
    LoadOpenGLExtension
  },
  {
    "TestActorLightingFlag",
    TestActorLightingFlag
  },
  {
    "TestAnimationScene",
    TestAnimationScene
  },
  {
    "TestBackfaceCulling",
    TestBackfaceCulling
  },
  {
    "TestBlurAndSobelPasses",
    TestBlurAndSobelPasses
  },
  {
    "TestDynamic2DLabelMapper",
    TestDynamic2DLabelMapper
  },
  {
    "TestFBO",
    TestFBO
  },
  {
    "TestFollowerPicking",
    TestFollowerPicking
  },
  {
    "TestGaussianBlurPass",
    TestGaussianBlurPass
  },
  {
    "TestGlyph3DMapper",
    TestGlyph3DMapper
  },
  {
    "TestGlyph3DMapperMasking",
    TestGlyph3DMapperMasking
  },
  {
    "TestGlyph3DMapperOrientationArray",
    TestGlyph3DMapperOrientationArray
  },
  {
    "TestGlyph3DMapperPicking",
    TestGlyph3DMapperPicking
  },
  {
    "TestGPUInfo",
    TestGPUInfo
  },
  {
    "TestGradientBackground",
    TestGradientBackground
  },
  {
    "TestHomogeneousTransformOfActor",
    TestHomogeneousTransformOfActor
  },
  {
    "TestImageResliceMapperAlpha",
    TestImageResliceMapperAlpha
  },
  {
    "TestImageResliceMapperBackground",
    TestImageResliceMapperBackground
  },
  {
    "TestImageResliceMapperBorder",
    TestImageResliceMapperBorder
  },
  {
    "TestImageResliceMapperInterpolation",
    TestImageResliceMapperInterpolation
  },
  {
    "TestImageResliceMapperOffAxis",
    TestImageResliceMapperOffAxis
  },
  {
    "TestImageResliceMapperOrient3D",
    TestImageResliceMapperOrient3D
  },
  {
    "TestImageResliceMapperSlab",
    TestImageResliceMapperSlab
  },
  {
    "TestImageSliceMapperAlpha",
    TestImageSliceMapperAlpha
  },
  {
    "TestImageSliceMapperBackground",
    TestImageSliceMapperBackground
  },
  {
    "TestImageSliceMapperBorder",
    TestImageSliceMapperBorder
  },
  {
    "TestImageSliceMapperOrient2D",
    TestImageSliceMapperOrient2D
  },
  {
    "TestImageSliceMapperOrient3D",
    TestImageSliceMapperOrient3D
  },
  {
    "TestImageSliceMapperInterpolation",
    TestImageSliceMapperInterpolation
  },
  {
    "TestImageStack",
    TestImageStack
  },
  {
    "TestInteractorTimers",
    TestInteractorTimers
  },
  {
    "TestLabelPlacer",
    TestLabelPlacer
  },
  {
    "TestLabelPlacer2D",
    TestLabelPlacer2D
  },
  {
    "TestLabelPlacerCoincidentPoints",
    TestLabelPlacerCoincidentPoints
  },
  {
    "TestLabelPlacementMapper2D",
    TestLabelPlacementMapper2D
  },
  {
    "TestLabelPlacementMapperCoincidentPoints",
    TestLabelPlacementMapperCoincidentPoints
  },
  {
    "TestLightActor",
    TestLightActor
  },
  {
    "TestLODActor",
    TestLODActor
  },
  {
    "TestManyActors",
    TestManyActors
  },
  {
    "TestOrderedTriangulator",
    TestOrderedTriangulator
  },
  {
    "TestOpacity",
    TestOpacity
  },
  {
    "TestOpenGLPolyDataMapper",
    TestOpenGLPolyDataMapper
  },
  {
    "TestOSConeCxx",
    TestOSConeCxx
  },
  {
    "TestPOVExporter",
    TestPOVExporter
  },
  {
    "TestResetCameraVerticalAspectRatio",
    TestResetCameraVerticalAspectRatio
  },
  {
    "TestResetCameraVerticalAspectRatioParallel",
    TestResetCameraVerticalAspectRatioParallel
  },
  {
    "TestSetImageOrientation",
    TestSetImageOrientation
  },
  {
    "TestSobelGradientMagnitudePass",
    TestSobelGradientMagnitudePass
  },
  {
    "TestShadowMapPass",
    TestShadowMapPass
  },
  {
    "TestTextActorAlphaBlending",
    TestTextActorAlphaBlending
  },
  {
    "TestTextActorDepthPeeling",
    TestTextActorDepthPeeling
  },
  {
    "TestTextActor3DAlphaBlending",
    TestTextActor3DAlphaBlending
  },
  {
    "TestTextActor3DDepthPeeling",
    TestTextActor3DDepthPeeling
  },
  {
    "TestTexturedBackground",
    TestTexturedBackground
  },
  {
    "TestTextureSize",
    TestTextureSize
  },
  {
    "TestTDx",
    TestTDx
  },
  {
    "TestTilingCxx",
    TestTilingCxx
  },
  {
    "TestTransformCoordinateUseDouble",
    TestTransformCoordinateUseDouble
  },
  {
    "TestTranslucentLUTAlphaBlending",
    TestTranslucentLUTAlphaBlending
  },
  {
    "TestTranslucentLUTDepthPeeling",
    TestTranslucentLUTDepthPeeling
  },
  {
    "TestTranslucentLUTDepthPeelingPass",
    TestTranslucentLUTDepthPeelingPass
  },
  {
    "TestTranslucentLUTTextureAlphaBlending",
    TestTranslucentLUTTextureAlphaBlending
  },
  {
    "TestTranslucentLUTTextureDepthPeeling",
    TestTranslucentLUTTextureDepthPeeling
  },
  {
    "TestGenericVertexAttributesGLSLCxx",
    TestGenericVertexAttributesGLSLCxx
  },
  {
    "TestGenericVertexAttributesGLSLAlphaBlending",
    TestGenericVertexAttributesGLSLAlphaBlending
  },
  {
    "TestGenericVertexAttributesGLSLDepthPeelingPass",
    TestGenericVertexAttributesGLSLDepthPeelingPass
  },

  { NULL, NULL } /* NOLINT */
};

static const int NumTests = CM_CAST(int,
  sizeof(cmakeGeneratedFunctionMapEntries) / sizeof(functionMapEntry)) - 1;

/* Allocate and create a lowercased copy of string
   (note that it has to be free'd manually) */
static char* lowercase(const char* string)
{
  char *new_string, *p;
  size_t stringSize;

  stringSize = CM_CAST(size_t, strlen(string) + 1);
  new_string = CM_CAST(char*, malloc(sizeof(char) * stringSize));

  if (new_string == NULL) { /* NOLINT */
    return NULL;            /* NOLINT */
  }
  strncpy(new_string, string, stringSize);
  for (p = new_string; *p != 0; ++p) {
    *p = CM_CAST(char, tolower(*p));
  }
  return new_string;
}

int main(int ac, char* av[])
{
  int i, testNum = 0, partial_match;
  char *arg, *test_name;
  int testToRun = -1;

  

  /* If no test name was given */
  /* process command line with user function.  */
  if (ac < 2) {
    /* Ask for a test.  */
    printf("Available tests:\n");
    for (i = 0; i < NumTests; ++i) {
      printf("%3d. %s\n", i, cmakeGeneratedFunctionMapEntries[i].name);
    }
    printf("To run a test, enter the test number: ");
    fflush(stdout);
    if (scanf("%d", &testNum) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
    if (testNum >= NumTests) {
      printf("%3d is an invalid test number.\n", testNum);
      return -1;
    }
    testToRun = testNum;
    ac--;
    av++;
  }
  partial_match = 0;
  arg = NULL; /* NOLINT */
  /* If partial match is requested.  */
  if (testToRun == -1 && ac > 1) {
    partial_match = (strcmp(av[1], "-R") == 0) ? 1 : 0;
  }
  if (partial_match != 0 && ac < 3) {
    printf("-R needs an additional parameter.\n");
    return -1;
  }
  if (testToRun == -1) {
    arg = lowercase(av[1 + partial_match]);
  }
  for (i = 0; i < NumTests && testToRun == -1; ++i) {
    test_name = lowercase(cmakeGeneratedFunctionMapEntries[i].name);
    if (partial_match != 0 && strstr(test_name, arg) != NULL) { /* NOLINT */
      testToRun = i;
      ac -= 2;
      av += 2;
    } else if (partial_match == 0 && strcmp(test_name, arg) == 0) {
      testToRun = i;
      ac--;
      av++;
    }
    free(test_name);
  }
  free(arg);
  if (testToRun != -1) {
    int result;

    vtkFloatingPointExceptions::Enable();

    try {
    if (testToRun < 0 || testToRun >= NumTests) {
      printf("testToRun was modified by TestDriver code to an invalid value: "
             "%3d.\n",
             testNum);
      return -1;
    }
    result = (*cmakeGeneratedFunctionMapEntries[testToRun].func)(ac, av);
    }
    catch(std::exception& e)
      {
      fprintf(stderr, "Test driver caught exception: [%s]\n", e.what());
      result = -1;
      }
    return result;
  }

  /* Nothing was run, display the test names.  */
  printf("Available tests:\n");
  for (i = 0; i < NumTests; ++i) {
    printf("%3d. %s\n", i, cmakeGeneratedFunctionMapEntries[i].name);
  }
  printf("Failed: %s is an invalid test name.\n", av[1]);

  return -1;
}
