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
	if ( argc != 25 ){
		cout << "USAGE: " << argv[0] << " tilt_angle axis_rotation axis_correction box_correction[2] normal[0] normal[1] normal[2] m[0-15] \n" << endl;
		return 0;
	}
#if 0
	Array< float, 3 > A(256,256,256);
	A = 0;
	//A( Range(49,50), 50, Range::all() ) = 1;
	//A(Range(127,128),Range(128),Range(127,128))=1;
	//A(Range(123,124),Range(120),Range(131,132))=1;
	A( 124, 127, 124 ) = 1;
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
	cout << "Tilt-angle, tilt axis=" << tilt_angle << "," << tilt_axis_angle << endl;

	double axis_correction = atof( argv[3] );
	double tilt_X_correction = atof( argv[4] );
	double tilt_Y_correction = atof( argv[5] );
	cout << "Shift corrections =[" << axis_correction << "," << tilt_X_correction << "," << tilt_Y_correction << "]" << endl;
	
	// spike normal
	double normal[] = { atof(argv[6]), atof(argv[7]), atof(argv[8])};
	cout << "Normal = " << normal[0] << "," << normal[1] << "," << normal[2] << endl;
	
	// refinement matrix
	double matrix[] = { atof(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atof(argv[14]), atof(argv[15]), atof(argv[16]), atof(argv[17]), atof(argv[18]), atof(argv[19]), atof(argv[20]), atof(argv[21]), atof(argv[22]), atof(argv[23]), atof(argv[24])};

	// refinement
	vtkMatrix4x4 * refinement = vtkMatrix4x4 :: New();
	refinement->DeepCopy( matrix );
	// cout << "Matrix = " << *refinement << endl;

	// net transformation matrix
	vtkTransform * t = vtkTransform :: New();

	// compensate for .box coordinate discretization
	t->Translate( tilt_X_correction, tilt_Y_correction, 0 );
	
	// convert to IMOD's 2D center of rotation
	t->Translate( .5, .5, 0 );
	
	// tilt axis angle rotation
	t->RotateZ( tilt_axis_angle );
	
	// undo center of rotation conversion
	t->Translate( -.5, -.5, 0 );

	// compute absolute orientation of tilt axis
	double nxtaxis [] = { 1, 0, 0, 1 };
	double nytaxis [] = { 0, 1, 0, 1 };
	double nztaxis [] = { 0, 0, 1, 1 };
	double xtaxis[4], ytaxis[4], ztaxis[4];
	
	vtkMatrix4x4 * inverse = vtkMatrix4x4 :: New();
	t->GetInverse( inverse );
	
	// get rid of translation component
	inverse->SetElement(0,3,0);
	inverse->SetElement(1,3,0);
	inverse->SetElement(2,3,0);
	
	// find vector in old XYZ coordinate system that gets maped into the new Y axis
	inverse->MultiplyPoint( nxtaxis, xtaxis );
	inverse->MultiplyPoint( nytaxis, ytaxis );
	inverse->MultiplyPoint( nztaxis, ztaxis );

	printf("xtaxis = [ %.5f, %.5f, %.5f ]\n",xtaxis[0],xtaxis[1],xtaxis[2]);
	printf("ytaxis = [ %.5f, %.5f, %.5f ]\n",ytaxis[0],ytaxis[1],ytaxis[2]);
	printf("ztaxis = [ %.5f, %.5f, %.5f ]\n",ztaxis[0],ztaxis[1],ztaxis[2]);
	
	// Correct for centering of 1D projection
	t->Translate( axis_correction * xtaxis[0], axis_correction * xtaxis[1], axis_correction * xtaxis[2] );

	// Convert to IMOD's tilt axis location
	double tx = .5 * fabs( xtaxis[0] );
	double ty = .5 * fabs( xtaxis[1] );
	double tz = .5;

	printf("\nTranslation to axis = [ %.5f %.5f %.5f ]\n", tx, ty, tz );

	t->Translate( tx, ty, tz );

	// apply tilt angle rotation
	t->RotateWXYZ( tilt_angle, ytaxis[0], ytaxis[1], ytaxis[2] );
	
	// undo tilt axis conversion
	t->Translate( -tx, -ty, -tz );
	
	// // Correct for centering of 1D projection
	// t->Translate( axis_correction, 0, 0 );
	
	// // translate to IMOD's tilt axis location
	// t->Translate( .5, 0, .5 );
	
	// // apply tilt angle rotation
	// t->RotateX( - tilt_angle );
	
	// // undo conversion
	// t->Translate( -.5, 0, -.5 );
	
	// cout << "Current t = "<< *t->GetMatrix() << endl;
	
	// apply spike euler angles
	vtkTransform * euler = vtkTransform::New();
	euler->RotateZ( - normal[2] );
	euler->RotateX( - normal[0] );
	euler->RotateZ( - normal[1] );
	
	// // tilt geometry
	// vtkTransform * tilt = vtkTransform::New();
	// tilt->RotateZ( tilt_axis_angle );
	// tilt->RotateY( - tilt_angle );
	// //tilt->RotateY( 180.0 );

	// t * euler * refinement
	vtkMatrix4x4 * product = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( t->GetMatrix(), euler->GetMatrix(), product );
	refinement->Delete();

	vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( product, refinement, final );
	
	// cout << "Final = " << *final << endl;
		
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
		double sign = final->GetElement(2,2) / fabs(final->GetElement(2,2));
		theta = vtkMath::Pi() * ( 1 - sign ) / 2.0;
		phi = 0.0;
		// THIS ORIGINAL CALCULATION WAS NOT CORRECT
		// psi = -	sign * atan2( final->GetElement(0,1) , final->GetElement(0,0) );
		float arg = final->GetElement(0,0) / cos( theta );
		psi = acos( arg );
	}

	// vtkTransform * t = vtkTransform::New();
	// t->RotateZ(phi* vtkMath::RadiansToDegrees());
	// t->RotateY(theta* vtkMath::RadiansToDegrees());
	// t->RotateZ(psi* vtkMath::RadiansToDegrees());
	// cout << "t = " << *t->GetMatrix() << endl;
 
	double frealign[3];
	// frealign[1] = acos( final->GetElement(2,2) ) * vtkMath::RadiansToDegrees();
	// frealign[0] = asin( final->GetElement(2,1) / sin( frealign[1] * vtkMath::DegreesToRadians() ) ) * vtkMath::RadiansToDegrees();
	// frealign[2] = asin( final->GetElement(1,2) / sin( frealign[1] * vtkMath::DegreesToRadians() ) ) * vtkMath::RadiansToDegrees();

	// // VTK's calculations
	// vtkTransform * tf = vtkTransform::New();
	// tf->SetMatrix(final);
	// float orient[3];
	// tf->GetOrientation(orient);
		
	frealign[0] = psi * vtkMath::RadiansToDegrees();
	frealign[1] = theta * vtkMath::RadiansToDegrees();
	frealign[2] = phi * vtkMath::RadiansToDegrees();
	
	// frealign does not use negative angles, so we add 360 to each negative angle
	frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
	frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
	frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];

	// Now project xyz shifts onto view plane for use in FREALIGN

	// retrieve shifts from final refinement
	double shifts [] = { final->GetElement(0,3), final->GetElement(1,3), final->GetElement(2,3) };
	
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
	
	double zaxis [] = { 0, 0, 1, 1 };
	double nzaxis[4];
	final->MultiplyPoint( zaxis, nzaxis );

	final->Delete();	

	double x3axis [] = { nxaxis[0], nxaxis[1], nxaxis[2] };
	double y3axis [] = { nyaxis[0], nyaxis[1], nyaxis[2] };
	double z3axis [] = { nzaxis[0], nzaxis[1], nzaxis[2] };
	
	// printf("nxaxis = [ %.5f, %.5f, %.5f ]\n",nxaxis[0],nxaxis[1],nxaxis[2]);
	// printf("nyaxis = [ %.5f, %.5f, %.5f ]\n",nyaxis[0],nyaxis[1],nyaxis[2]);
	// printf("nzaxis = [ %.5f, %.5f, %.5f ]\n",nzaxis[0],nzaxis[1],nzaxis[2]);

	double sx = vtkMath::Dot( shifts, x3axis );
	double sy = vtkMath::Dot( shifts, y3axis );
	double sz = vtkMath::Dot( shifts, z3axis );

	cout << "[Sx,Sy,Sz]=[" << sx << "," << sy << "," << sz << "]" << endl;

	/////////////////////////////////////////////////////////////////////////////
	// WARNING! BEFORE MEMORIAL DAY WEEKEND I CHANGED THE ORDER OF PHI AND PSI //
	/////////////////////////////////////////////////////////////////////////////
	printf("\nParameters for FREALIGN [PHI,THETA,PSI,SX,SY]= [ %.5f %.5f %.5f %.5f %.5f ]\n", frealign[2], frealign[1], frealign[0], sx, sy );
	
	return 0;
}
