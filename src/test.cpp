
#define NBF_PARALLEL_IMPLEMENTATION_MPI 1
#define NBF_AVERAGE_IN_RECIPROCAL_SPACE
#define BZ_GENERATE_GLOBAL_INSTANCES

#include "mpi.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <blitz/bzdebug.h>        // Test suite globals

using namespace blitz;

#include <string.h>

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfUnsupervisedLoopClustering.h>

#include <nbfCylindricalDomain3.h>

#include <mxml.h>

extern "C" {
#include <svdlib.h>
}

#include <nbfTimer.h>

#define PIXEL float



#include <iostream>
#include <stdio.h>
//using namespace std;


//static nbfFourierFilter< PIXEL, 3 > g_fourierFilter;

class Tester {
	public:
		Tester();
	private:
		vtkImageFourierCenter * p;
};

Tester::Tester() {
	vtkImageFourierCenter::New();
}

static Tester tester;


int main( int argc, char ** argv )
{
	// TEMP
	cout << "A" << endl;
	
	//nbfFourierFilter< PIXEL, 3 > fourierFilter;
	//void * p = &g_fourierFilter;
	//void * p = &tester;
	//cout << p << endl;

	vtkImageFourierCenter * thing = vtkImageFourierCenter::New();
	cout << thing << endl;

	//nbfWedgedAverageImage3D< PIXEL > foo;
	//nbfWedgedSubImage3D< PIXEL > foo;

	// TEMP
	cout << "B" << endl;

	return 0;
}
