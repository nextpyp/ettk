#include "mpi.h"
#define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

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
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkPoints.h>
#include <vtkPolyVertex.h>
#include <vtkProbeFilter.h>
#include <vtkContourFilter.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMrcWriter.h>

#include <em/nbfFourierImageMetricCore.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfExtractPointsAndNormals3D.h>
#include <em/nbfTemplateSearchEM.h>
#include <em/nbfCorrelationImageMetric.h>

#include <em/nbfCutSubVolumes.h>


#define PIXEL double

int main( int argc, char ** argv )
{
	if ( argc != 22 ){
		cout << "USAGE: " << argv[0] << " tilt_angle axis_rotation normal[0] normal[1] normal[2] m[0-15] \n" << endl;
		return 0;
	}
#if 0
	Array< float, 3 > A(100,100,90);
	A = 0;
	A( Range(49,50), 50, Range::all() ) = 1;
	//A( 50, 50, 50 ) = 1;
	vtkImageData * data = vtkImageData::New();
	//nbfMrcReader reader;
	//reader.setFileName( argv[1] );
	//reader.read(data);
	nbfVTKInterface::blitzToVtk(A,data);
	nbfMrcWriter w;
	w.setFileName("phantom.mrc");
	w.write(data);
	return 0;
#endif	
	
	// tilt angle
	double tilt_angle = atof( argv[1] );

	// tilt axis rotation
	double tilt_axis_angle = atof( argv[2] );

	//cout << "Tilt-angle, tilt axis=" << tilt_angle << "," << tilt_axis_angle << endl;
	
	// tilt geometry
	vtkTransform * tilt = vtkTransform::New();
	tilt->RotateZ( tilt_axis_angle );
	tilt->RotateY( - tilt_angle );
	//tilt->RotateY( 180.0 );

	// IMOD-FREALIGN geometry adaptation
	vtkTransform * t = vtkTransform::New();
	//t->Translate( .5, 0, .5 );

	// spike normal
	double normal[] = { atof(argv[3]), atof(argv[4]), atof(argv[5])};
	cout << "Normal = " << normal[0] << "," << normal[1] << "," << normal[2] << endl;
	
	vtkTransform * euler = vtkTransform::New();
	euler->RotateZ( normal[2] );
	euler->RotateX( normal[0] );
	euler->RotateZ( normal[1] );

	// t * euler
	vtkMatrix4x4 * teuler = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( t->GetMatrix(), euler->GetMatrix(), teuler );
	euler->Delete();
	t->Delete();
	
	// refinement matrix
	double matrix[] = { atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atof(argv[14]), atof(argv[15]), atof(argv[16]), atof(argv[17]), atof(argv[18]), atof(argv[19]), atof(argv[20]), atof(argv[21])};

	// refinement
	vtkMatrix4x4 * refinement = vtkMatrix4x4 :: New();
	refinement->DeepCopy( matrix );

	cout << "Matrix = " << *refinement << endl;

	// teuler * refinement
	vtkMatrix4x4 * product = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( teuler, refinement, product );
	refinement->Delete();

	vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( tilt->GetMatrix(), product, final );
	product->Delete();		
	tilt->Delete();

	// // for some reason the spike ends upside down. This transformation will fix the orientation of the Z axis.
	// vtkTransform * c = vtkTransform :: New();
	// //c->RotateX(-180);
	// //c->RotateZ(-180);
	// //c->RotateY(180);
	// c->Translate( -.5, 0, .5 );
	// vtkMatrix4x4::Multiply4x4( final, c->GetMatrix(), product );
	// final->DeepCopy( product );
	// product->Delete();
	// c->Delete();
	
	cout << "Final = " << *final << endl;
		
	// my calculations
	double psi, theta, phi;
	
	if ( fabs(final->GetElement(2,2)) < 1 - numeric_limits< float > :: min() ){
		theta = acos( final->GetElement(2,2) );
		double num = final->GetElement(2,1) / sin( theta );
		double den = - final->GetElement(2,0) / sin( theta );
		psi = atan2( num, den );
	
		num = final->GetElement(1,2) / sin( theta );
		den = final->GetElement(0,2) / sin( theta );
		phi = atan2( num, den );
	} else {
		cout << "Using alternate" << endl;
		//cout << "Final = " << *final << endl;
		double sign = final->GetElement(2,2) / fabs(final->GetElement(2,2));
		theta = vtkMath::Pi() * ( 1 - sign ) / 2.0;
		phi = 0.0;
		// THIS ORIGINAL CALCULATION WAS NOT CORRECT
		// psi = -	sign * atan2( final->GetElement(0,1) , final->GetElement(0,0) );
		float arg = final->GetElement(0,0) / cos( theta );
		psi = acos( arg );
	}

	// cout << "[" << theta1 << "," << psi1 << "," << phi1 << "]" << endl;
	// cout << "[" << theta2 << "," << psi2 << "," << phi2 << "]" << endl;
	
	// vtkTransform * t1 = vtkTransform::New();
	// t1->RotateZ(phi* vtkMath::RadiansToDegrees());
	// t1->RotateY(theta* vtkMath::RadiansToDegrees());
	// t1->RotateZ(psi* vtkMath::RadiansToDegrees());
	// cout << *t1->GetMatrix() << endl;
 
    // vtkTransform * t1 = vtkTransform::New();
	// t1->RotateZ(phi* vtkMath::RadiansToDegrees());
	// t1->RotateY(theta* vtkMath::RadiansToDegrees());
	// t1->RotateZ(psi* vtkMath::RadiansToDegrees());
	// cout << *t1->GetMatrix() << endl;

	// vtkTransform * t2 = vtkTransform::New();
	// t2->RotateZ(phi2* vtkMath::RadiansToDegrees());
	// t2->RotateY(theta2* vtkMath::RadiansToDegrees());
	// t2->RotateZ(psi2* vtkMath::RadiansToDegrees());
	// cout << *t2->GetMatrix() << endl;
	
	double frealign[3];
	// frealign[1] = acos( final->GetElement(2,2) ) * vtkMath::RadiansToDegrees();
	// frealign[0] = asin( final->GetElement(2,1) / sin( frealign[1] * vtkMath::DegreesToRadians() ) ) * vtkMath::RadiansToDegrees();
	// frealign[2] = asin( final->GetElement(1,2) / sin( frealign[1] * vtkMath::DegreesToRadians() ) ) * vtkMath::RadiansToDegrees();

	// // VTK's calculations
	// vtkTransform * tf = vtkTransform::New();
	// tf->SetMatrix(final);
	// float orient[3];
	// tf->GetOrientation(orient);
		
	// using VTK's calculations
	frealign[0] = vtkMath::DegreesFromRadians(psi);
	frealign[1] = vtkMath::DegreesFromRadians(theta);
	frealign[2] = vtkMath::DegreesFromRadians(phi);
	
	// frealign does not use negative angles, so we add 360 to each negative angle
	frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
	frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
	frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];

	// Now project xyz shifts onto view plane for use in FREALIGN

	// retrieve shifts from final refinement
	// double shifts [] = { atof(argv[9]), atof(argv[13]), atof(argv[17]) };
	//double shifts [] = { matrix[3], matrix[7], matrix[11] };
	double shifts [] = { product->GetElement(0,3), product->GetElement(1,3), product->GetElement(2,3) };
	
	cout << "S=[" << shifts[0] << "," << shifts[1] << "," << shifts[2] << "]" << endl;

	// set translation component to 0
	final->SetElement(0,3,0);
	final->SetElement(1,3,0);
	final->SetElement(2,3,0);
	
	double xaxis [] = { 1, 0, 0, 1 };
	double nxaxis[4];
	final->MultiplyPoint( xaxis, nxaxis );
		
	double yaxis [] = { 0, 1, 0, 1 };
	double nyaxis[4];
	final->MultiplyPoint( yaxis, nyaxis );
	
	double x3axis [] = { nxaxis[0], nxaxis[1], nxaxis[2] };
	double y3axis [] = { nyaxis[0], nyaxis[1], nyaxis[2] };
	
	double sx = vtkMath::Dot( shifts, x3axis );
	double sy = vtkMath::Dot( shifts, y3axis );

	// cout << x3axis[0] << "," << x3axis[1] << "," << x3axis[2] << endl;
	// cout << y3axis[0] << "," << y3axis[1] << "," << y3axis[2] << endl;
	
	// final->MultiplyPoint( yaxis, nyaxis );
	// cout << nyaxis[0] << "," << nyaxis[1] << "," << nyaxis[2] << endl;
	// t1->GetMatrix()->MultiplyPoint( yaxis, nyaxis );
	// cout << nyaxis[0] << "," << nyaxis[1] << "," << nyaxis[2] << endl;
	// t2->GetMatrix()->MultiplyPoint( yaxis, nyaxis );
	// cout << nyaxis[0] << "," << nyaxis[1] << "," << nyaxis[2] << endl;

	printf("\nParameters for FREALIGN [PHI,THETA,PSI,SX,SY]= [ %.5f %.5f %.5f %.5f %.5f ]\n", frealign[0], frealign[1], frealign[2], sx, sy );
	
	//cout << "Parameters for FREALIGN [PHI,THETA,PSI,SX,SY]= [ " << frealign[0] << " " << frealign[1] << " " << frealign[2] << " " << sx << " " << sy << " ]" << endl;
	
	final->Delete();
	
	return 0;
}
