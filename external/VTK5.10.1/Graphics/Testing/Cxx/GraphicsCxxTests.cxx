#include <ctype.h>  /* NOLINT */
#include <stdio.h>  /* NOLINT */
#include <stdlib.h> /* NOLINT */
#include <string.h> /* NOLINT */

#if defined(_MSC_VER)
#pragma warning(disable : 4996) /* deprecation */
#endif

#include "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Rendering/vtkTestingObjectFactory.h"


/* Forward declare test functions. */
int TestCenterOfMass(int, char*[]);
int Mace(int, char*[]);
int expCos(int, char*[]);
int BoxClipTriangulate(int, char*[]);
int CellLocator(int, char*[]);
int PointLocator(int, char*[]);
int FrustumClip(int, char*[]);
int RGrid(int, char*[]);
int TestAppendSelection(int, char*[]);
int TestAssignAttribute(int, char*[]);
int TestBareScalarsToColors(int, char*[]);
int TestBSPTree(int, char*[]);
int TestCellDataToPointData(int, char*[]);
int TestDensifyPolyData(int, char*[]);
int TestClipHyperOctree(int, char*[]);
int TestConvertSelection(int, char*[]);
int TestDelaunay2D(int, char*[]);
int TestGlyph3D(int, char*[]);
int TestExtraction(int, char*[]);
int TestExtractSelection(int, char*[]);
int TestExtractSurfaceNonLinearSubdivision(int, char*[]);
int TestHyperOctreeContourFilter(int, char*[]);
int TestHyperOctreeCutter(int, char*[]);
int TestHyperOctreeDual(int, char*[]);
int TestHyperOctreeSurfaceFilter(int, char*[]);
int TestHyperOctreeToUniformGrid(int, char*[]);
int TestImageDataToPointSet(int, char*[]);
int TestLineSource(int, char*[]);
int TestMapVectorsAsRGBColors(int, char*[]);
int TestMapVectorsToColors(int, char*[]);
int TestNamedComponents(int, char*[]);
int TestMeanValueCoordinatesInterpolation1(int, char*[]);
int TestMeanValueCoordinatesInterpolation2(int, char*[]);
int TestPolyDataPointSampler(int, char*[]);
int TestPolyhedron0(int, char*[]);
int TestPolyhedron1(int, char*[]);
int TestQuadRotationalExtrusion(int, char*[]);
int TestRectilinearGridToPointSet(int, char*[]);
int TestReflectionFilter(int, char*[]);
int TestRotationalExtrusion(int, char*[]);
int TestSelectEnclosedPoints(int, char*[]);
int TestTessellatedBoxSource(int, char*[]);
int TestTessellator(int, char*[]);
int TestUncertaintyTubeFilter(int, char*[]);
int TestDecimatePolylineFilter(int, char*[]);
int TestImplicitPolyDataDistance(int, char*[]);


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
    "TestCenterOfMass",
    TestCenterOfMass
  },
  {
    "Mace",
    Mace
  },
  {
    "expCos",
    expCos
  },
  {
    "BoxClipTriangulate",
    BoxClipTriangulate
  },
  {
    "CellLocator",
    CellLocator
  },
  {
    "PointLocator",
    PointLocator
  },
  {
    "FrustumClip",
    FrustumClip
  },
  {
    "RGrid",
    RGrid
  },
  {
    "TestAppendSelection",
    TestAppendSelection
  },
  {
    "TestAssignAttribute",
    TestAssignAttribute
  },
  {
    "TestBareScalarsToColors",
    TestBareScalarsToColors
  },
  {
    "TestBSPTree",
    TestBSPTree
  },
  {
    "TestCellDataToPointData",
    TestCellDataToPointData
  },
  {
    "TestDensifyPolyData",
    TestDensifyPolyData
  },
  {
    "TestClipHyperOctree",
    TestClipHyperOctree
  },
  {
    "TestConvertSelection",
    TestConvertSelection
  },
  {
    "TestDelaunay2D",
    TestDelaunay2D
  },
  {
    "TestGlyph3D",
    TestGlyph3D
  },
  {
    "TestExtraction",
    TestExtraction
  },
  {
    "TestExtractSelection",
    TestExtractSelection
  },
  {
    "TestExtractSurfaceNonLinearSubdivision",
    TestExtractSurfaceNonLinearSubdivision
  },
  {
    "TestHyperOctreeContourFilter",
    TestHyperOctreeContourFilter
  },
  {
    "TestHyperOctreeCutter",
    TestHyperOctreeCutter
  },
  {
    "TestHyperOctreeDual",
    TestHyperOctreeDual
  },
  {
    "TestHyperOctreeSurfaceFilter",
    TestHyperOctreeSurfaceFilter
  },
  {
    "TestHyperOctreeToUniformGrid",
    TestHyperOctreeToUniformGrid
  },
  {
    "TestImageDataToPointSet",
    TestImageDataToPointSet
  },
  {
    "TestLineSource",
    TestLineSource
  },
  {
    "TestMapVectorsAsRGBColors",
    TestMapVectorsAsRGBColors
  },
  {
    "TestMapVectorsToColors",
    TestMapVectorsToColors
  },
  {
    "TestNamedComponents",
    TestNamedComponents
  },
  {
    "TestMeanValueCoordinatesInterpolation1",
    TestMeanValueCoordinatesInterpolation1
  },
  {
    "TestMeanValueCoordinatesInterpolation2",
    TestMeanValueCoordinatesInterpolation2
  },
  {
    "TestPolyDataPointSampler",
    TestPolyDataPointSampler
  },
  {
    "TestPolyhedron0",
    TestPolyhedron0
  },
  {
    "TestPolyhedron1",
    TestPolyhedron1
  },
  {
    "TestQuadRotationalExtrusion",
    TestQuadRotationalExtrusion
  },
  {
    "TestRectilinearGridToPointSet",
    TestRectilinearGridToPointSet
  },
  {
    "TestReflectionFilter",
    TestReflectionFilter
  },
  {
    "TestRotationalExtrusion",
    TestRotationalExtrusion
  },
  {
    "TestSelectEnclosedPoints",
    TestSelectEnclosedPoints
  },
  {
    "TestTessellatedBoxSource",
    TestTessellatedBoxSource
  },
  {
    "TestTessellator",
    TestTessellator
  },
  {
    "TestUncertaintyTubeFilter",
    TestUncertaintyTubeFilter
  },
  {
    "TestDecimatePolylineFilter",
    TestDecimatePolylineFilter
  },
  {
    "TestImplicitPolyDataDistance",
    TestImplicitPolyDataDistance
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

    // Set defaults
    vtkTestingInteractor::ValidBaseline =
      std::string("VTK_DATA_ROOT-NOTFOUND") +
      std::string("/Baseline/") +
      std::string("Graphics/") +
      std::string(cmakeGeneratedFunctionMapEntries[testToRun].name) +
      std::string(".png");
    vtkTestingInteractor::TempDirectory =
      std::string("/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Testing/Temporary");
    vtkTestingInteractor::DataDirectory =
      std::string("VTK_DATA_ROOT-NOTFOUND");

    int interactive = 0;
    for (int ii = 0; ii < ac; ++ii)
      {
      if ( strcmp(av[ii],"-I") == 0)
        {
        interactive = 1;
        continue;
        }
      if ( strcmp(av[ii],"-V") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::ValidBaseline = std::string(av[ii+1]);
        ++ii;
        continue;
        }
      if ( strcmp(av[ii],"-T") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::TempDirectory = std::string(av[ii+1]);
        ++ii;
        continue;
        }
      if ( strcmp(av[ii],"-D") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::DataDirectory = std::string(av[ii+1]);
        ++ii;
        continue;
        }
      if ( strcmp(av[ii],"-E") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::ErrorThreshold =
          static_cast<double>(atof(av[ii+1]));
        ++ii;
        continue;
        }
      }
    vtkSmartPointer<vtkTestingObjectFactory> factory = vtkSmartPointer<vtkTestingObjectFactory>::New();
    if (!interactive)
      {
      vtkObjectFactory::RegisterFactory(factory);
      }

    if (testToRun < 0 || testToRun >= NumTests) {
      printf("testToRun was modified by TestDriver code to an invalid value: "
             "%3d.\n",
             testNum);
      return -1;
    }
    result = (*cmakeGeneratedFunctionMapEntries[testToRun].func)(ac, av);
    
   if (!interactive)
     {
     if (vtkTestingInteractor::TestReturnStatus != -1)
        {
        if( vtkTestingInteractor::TestReturnStatus != vtkTesting::PASSED)
          {
          result = EXIT_FAILURE;
          }
        else
          {
          result = EXIT_SUCCESS;
          }
        }
      vtkObjectFactory::UnRegisterFactory(factory);
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
