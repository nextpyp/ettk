#include "mpi.h"
#define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES


#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <limits>

#define PIXEL double

using namespace std;

int main( int argc, char ** argv )
{
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	if ( argc != 25 ){
		cout << "USAGE: " << argv[0] << " tilt_angle axis_rotation box_correction[2] normal[0] normal[1] normal[2] m[0-15] cutOffset \n" << endl;
		return 0;
	}
	// tilt angle
	double tilt_angle = atof( argv[1] );

	// tilt axis rotation
	double tilt_axis_angle = atof( argv[2] );

	double tilt_X_correction = atof( argv[3] );
	double tilt_Y_correction = atof( argv[4] );
	
	// spike normal
	double normal[] = { atof(argv[5]), atof(argv[6]), atof(argv[7])};

	cout << "Tilt axis angle = " << tilt_axis_angle << endl;
	cout << "Tilt angle = " << tilt_angle << endl;
	cout << "Shift corrections =[" << tilt_X_correction << "," << tilt_Y_correction << "]" << endl;
	cout << "Normal = " << normal[0] << "," << normal[1] << "," << normal[2] << endl;
	
	// 3DAVG refinement
	// This transformation matrix is different from the matrix[0-15] in the 3DAVG refinement file. 
	// The ordering of the rotations has to be changed to match the new geometry. 
	// 	- The ordering in 3DAVG is: RotZ1 * RotX * RotZ2.
	//	- The ordering here is: RotZ2 * RotX * RotZ1
	// The new transformation matrix has the following expression:

	double matrix[] = { atof(argv[8]), -atof(argv[12]), atof(argv[16]), atof(argv[11]), 
	                   -atof(argv[9]), atof(argv[13]),-atof(argv[17]), -atof(argv[15]),
						atof(argv[10]),-atof(argv[14]), atof(argv[18]), atof(argv[19]),
						atof(argv[20]), atof(argv[21]), atof(argv[22]), atof(argv[23])};

	double rotationmatrix[] = { atof(argv[8]), -atof(argv[12]), atof(argv[16]), 0, 
							   -atof(argv[9]), atof(argv[13]),-atof(argv[17]), 0,
						        atof(argv[10]),-atof(argv[14]), atof(argv[18]), 0,
						        atof(argv[20]), atof(argv[21]), atof(argv[22]), atof(argv[23])};

	double cutOffset = atof( argv[24] );
	
	// 3DAVG refinement transformation
	vtkMatrix4x4 * refinement = vtkMatrix4x4 :: New();
	refinement->DeepCopy( matrix );

	// 3DAVG refinement transformation (rotation component only)
	vtkMatrix4x4 * refinementRotation = vtkMatrix4x4 :: New();
	refinementRotation->DeepCopy( rotationmatrix );
	
	// transformation matrix
	vtkTransform * t = vtkTransform :: New();

	// tilt axis angle rotation
	t->RotateZ( -tilt_axis_angle );

	// this is the mismatch between the reconstructions (need to double check z coordinate to make sure it's -.5 and not .5)
	t->Translate( .5, 0, -.5 );

	// apply tilt angle rotation
	t->RotateY( tilt_angle );

	// The remaining transformations are the spike normal and the 3DAVG refinement transformation.
	// Spike normals are a pure rotation R1. 3DAVG refinement is a full rotation and translation matrix F = R2 * T2
	// If the two transformations are composed, the net rotation is then: R = R1 * R2 and the transaltion component is T2.
	// The origin of the net rotation R is different from the image rotation, so to account for this we express
	// as R * T3, where T3 corrects for the difference in the origin of the rotation.
	
	// apply spike euler angles (pure rotation R1)
	t->RotateZ( normal[2] );
	t->RotateX( normal[0] );
	t->RotateZ( normal[1] );

	vtkMatrix4x4 * local = vtkMatrix4x4 :: New();
	
	// apply 3DAVG refinement transformation (rotation only)
	vtkMatrix4x4 :: Multiply4x4( t->GetMatrix(), refinementRotation, local );
	t->SetMatrix( local );

	// Difference vector between rotation origins C1=[51,50,50] and C2=[51,51,50]
	double diff [] = { 0, 1, 0, 1 };

	vtkTransform * r = vtkTransform :: New();
	r->RotateZ( normal[2] );
	r->RotateX( normal[0] );
	r->RotateZ( normal[1] );
	// apply 3DAVG refinement transformation (rotation only)
	vtkMatrix4x4::Multiply4x4( r->GetMatrix(), refinementRotation, local );
	r->SetMatrix( local );
	
	// Compute: t = Rot(C1) wrt C2 - C1 = Rot * ( C1 - C2 ) - ( C1 - C2 )
	double tdiff[4];	
	r->MultiplyPoint( diff, tdiff );
	printf("tdiff = [ %.5f, %.5f, %.5f ]\n",tdiff[0],tdiff[1],tdiff[2]);
	tdiff[0] -= diff[0];
	tdiff[1] -= diff[1];
	tdiff[2] -= diff[2];
	
	printf("tdiff = [ %.5f, %.5f, %.5f ]\n",tdiff[0],tdiff[1],tdiff[2]);

	// compute post-multiplying translation component from refinement matrix, T2
	refinementRotation->Invert();
	vtkMatrix4x4 * f = vtkMatrix4x4 :: New();
	vtkMatrix4x4::Multiply4x4( refinementRotation, refinement, f );

	printf("f translations = [ %.5f, %.5f, %.5f ]\n",f->GetElement(0,3),f->GetElement(1,3),f->GetElement(2,3));

	t->Translate( -tdiff[0] + f->GetElement(0,3), -tdiff[1] - f->GetElement(1,3), -tdiff[2] + f->GetElement(2,3) );
	
	// apply the corrections due to .box discretization (on the transformed coordinate system)
	t->Translate( -tilt_X_correction, tilt_Y_correction, 0 );

	vtkMatrix4x4 * final = t->GetMatrix();
	
	cout << "final = " <<  *final << endl;

	// extract euler angles from transformation matrix
	double psi, theta, phi;
	
	if ( fabs(final->GetElement(2,2)) < 1 - numeric_limits< float > :: min() ){
		theta = acos( final->GetElement(2,2) );
		double num = final->GetElement(1,2) / sin( theta );
		double den = - final->GetElement(0,2) / sin( theta );
		psi = atan2( num, den );
	
		num = final->GetElement(2,1) / sin( theta );
		den = final->GetElement(2,0) / sin( theta );
		phi = atan2( num, den );
	} else {
		double sign = final->GetElement(2,2) / fabs(final->GetElement(2,2));
		theta = vtkMath::Pi() * ( 1 - sign ) / 2.0;
		phi = 0.0;
		// THIS ORIGINAL CALCULATION WAS NOT CORRECT
		// psi = -	sign * atan2( final->GetElement(0,1) , final->GetElement(0,0) );
		float arg = final->GetElement(0,0) / cos( theta );
		psi = acos( arg );
	}

	double frealign[3];
		
	frealign[0] = psi * vtkMath::RadiansToDegrees();
	frealign[1] = theta * vtkMath::RadiansToDegrees();
	frealign[2] = phi * vtkMath::RadiansToDegrees();
	
	// frealign does not use negative angles, so we add 360 to each negative angle
	frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
	frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
	frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];

	// Now project xyz shifts onto view plane for use in FREALIGN

	// compute inverse of FREALIGN rotation transformation
	vtkTransform * rt = vtkTransform::New();
	rt->RotateZ( phi * vtkMath::RadiansToDegrees() );
	rt->RotateY( theta * vtkMath::RadiansToDegrees() );
	rt->RotateZ( psi * vtkMath::RadiansToDegrees() );

	vtkMatrix4x4 * nt = vtkMatrix4x4::New();
	vtkMatrix4x4::Multiply4x4( rt->GetMatrix(), final, nt );

	// printf( "final = [ %f, %f, %f ]\n", final->GetElement(0,3),final->GetElement(1,3),final->GetElement(2,3) );
	// printf( "nt = [ %f, %f, %f ]\n", nt->GetElement(0,3),nt->GetElement(1,3),nt->GetElement(2,3) );

	// vtkMatrix4x4::Multiply4x4( final, l, nt );
	// printf( "inv nt = [ %f, %f, %f ]\n", nt->GetElement(0,3),nt->GetElement(1,3),nt->GetElement(2,3) );
	// cout << *nt << endl;
	
	double sx = nt->GetElement(0,3);
	double sy = nt->GetElement(1,3);
	double sz = nt->GetElement(2,3);
	
	printf("\nParameters for FREALIGN [PHI,THETA,PSI,SX,SY]= [ %.5f %.5f %.5f %.5f %.5f ]\n", frealign[2], frealign[1], frealign[0], sx, sy );
	
	return 0;
}
