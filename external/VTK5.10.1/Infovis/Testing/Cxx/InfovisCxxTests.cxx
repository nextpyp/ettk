#include <ctype.h>  /* NOLINT */
#include <stdio.h>  /* NOLINT */
#include <stdlib.h> /* NOLINT */
#include <string.h> /* NOLINT */

#if defined(_MSC_VER)
#pragma warning(disable : 4996) /* deprecation */
#endif

#include "vtkTestDriver.h"


/* Forward declare test functions. */
int TestArcEdges(int, char*[]);
int TestAssignCoordinates(int, char*[]);
int TestCirclePackLayoutStrategy(int, char*[]);
int TestCollapseVerticesByArray(int, char*[]);
int TestContingencyStatistics(int, char*[]);
int TestConvexHull2D(int, char*[]);
int TestCorrelativeStatistics(int, char*[]);
int TestCosmicTreeLayoutStrategy(int, char*[]);
int TestDataObjectToTable(int, char*[]);
int TestDescriptiveStatistics(int, char*[]);
int TestExtractSelectedGraph(int, char*[]);
int TestGraph(int, char*[]);
int TestGraphAttributes(int, char*[]);
int TestGraphAlgorithms(int, char*[]);
int TestGraphHierarchicalBundle(int, char*[]);
int TestGraphLayoutStrategy(int, char*[]);
int TestGroupLeafVertices(int, char*[]);
int TestKCoreDecomposition(int, char*[]);
int TestKdTreeBoxSelection(int, char*[]);
int TestKMeansStatistics(int, char*[]);
int TestMergeGraphs(int, char*[]);
int TestMultiCorrelativeStatistics(int, char*[]);
int TestOrderStatistics(int, char*[]);
int TestPassArrays(int, char*[]);
int TestPassThrough(int, char*[]);
int TestPCAStatistics(int, char*[]);
int TestPruneTreeFilter(int, char*[]);
int TestRandomGraphSource(int, char*[]);
int TestRemoveIsolatedVertices(int, char*[]);
int TestSimple3DCirclesStrategy(int, char*[]);
int TestStreamGraph(int, char*[]);
int TestStringToCategory(int, char*[]);
int TestTable(int, char*[]);
int TestTableSplitColumnComponents(int, char*[]);
int TestTreeMapLayoutStrategy(int, char*[]);
int TestThresholdTable(int, char*[]);
int TestVariantArray(int, char*[]);
int TestVariant(int, char*[]);
int TestSQLGraphReader(int, char*[]);
int TestTimePoint(int, char*[]);


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
    "TestArcEdges",
    TestArcEdges
  },
  {
    "TestAssignCoordinates",
    TestAssignCoordinates
  },
  {
    "TestCirclePackLayoutStrategy",
    TestCirclePackLayoutStrategy
  },
  {
    "TestCollapseVerticesByArray",
    TestCollapseVerticesByArray
  },
  {
    "TestContingencyStatistics",
    TestContingencyStatistics
  },
  {
    "TestConvexHull2D",
    TestConvexHull2D
  },
  {
    "TestCorrelativeStatistics",
    TestCorrelativeStatistics
  },
  {
    "TestCosmicTreeLayoutStrategy",
    TestCosmicTreeLayoutStrategy
  },
  {
    "TestDataObjectToTable",
    TestDataObjectToTable
  },
  {
    "TestDescriptiveStatistics",
    TestDescriptiveStatistics
  },
  {
    "TestExtractSelectedGraph",
    TestExtractSelectedGraph
  },
  {
    "TestGraph",
    TestGraph
  },
  {
    "TestGraphAttributes",
    TestGraphAttributes
  },
  {
    "TestGraphAlgorithms",
    TestGraphAlgorithms
  },
  {
    "TestGraphHierarchicalBundle",
    TestGraphHierarchicalBundle
  },
  {
    "TestGraphLayoutStrategy",
    TestGraphLayoutStrategy
  },
  {
    "TestGroupLeafVertices",
    TestGroupLeafVertices
  },
  {
    "TestKCoreDecomposition",
    TestKCoreDecomposition
  },
  {
    "TestKdTreeBoxSelection",
    TestKdTreeBoxSelection
  },
  {
    "TestKMeansStatistics",
    TestKMeansStatistics
  },
  {
    "TestMergeGraphs",
    TestMergeGraphs
  },
  {
    "TestMultiCorrelativeStatistics",
    TestMultiCorrelativeStatistics
  },
  {
    "TestOrderStatistics",
    TestOrderStatistics
  },
  {
    "TestPassArrays",
    TestPassArrays
  },
  {
    "TestPassThrough",
    TestPassThrough
  },
  {
    "TestPCAStatistics",
    TestPCAStatistics
  },
  {
    "TestPruneTreeFilter",
    TestPruneTreeFilter
  },
  {
    "TestRandomGraphSource",
    TestRandomGraphSource
  },
  {
    "TestRemoveIsolatedVertices",
    TestRemoveIsolatedVertices
  },
  {
    "TestSimple3DCirclesStrategy",
    TestSimple3DCirclesStrategy
  },
  {
    "TestStreamGraph",
    TestStreamGraph
  },
  {
    "TestStringToCategory",
    TestStringToCategory
  },
  {
    "TestTable",
    TestTable
  },
  {
    "TestTableSplitColumnComponents",
    TestTableSplitColumnComponents
  },
  {
    "TestTreeMapLayoutStrategy",
    TestTreeMapLayoutStrategy
  },
  {
    "TestThresholdTable",
    TestThresholdTable
  },
  {
    "TestVariantArray",
    TestVariantArray
  },
  {
    "TestVariant",
    TestVariant
  },
  {
    "TestSQLGraphReader",
    TestSQLGraphReader
  },
  {
    "TestTimePoint",
    TestTimePoint
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
