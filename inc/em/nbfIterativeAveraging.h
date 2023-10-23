#pragma once

using namespace blitz;

#include <vector>
#include <em/nbfImageMetric.h>

/** Interface for Iterative Averaging Methods.
*/
template< class Pixel >
class nbfIterativeAveraging
{
public:

	nbfIterativeAveraging();

	~nbfIterativeAveraging(){};

	/// Set input list of images.
	void setInput( vector< vtkImageData * > & in ){ this->input = in; } 

	/// Set metric to use in computations.
	void setMetric( nbfImageMetric< Pixel, 3 > * m ){ this->metric = m; }

	/// Run averaging algorithm for a number of iterations.
	void execute( int, vtkImageData * );

protected:

	vector< vtkImageData * > input;
	nbfImageMetric< Pixel, 3 > * metric;
};

template< class Pixel >
nbfIterativeAveraging< Pixel > :: nbfIterativeAveraging()
{
	this->metric = NULL;
}

template< class Pixel >
void nbfIterativeAveraging< Pixel > :: execute( int numberOfIterations, vtkImageData * output )
{
	// average image
	vtkImageData * average = vtkImageData::New();

	// initialize to first image
	average->ShallowCopy( this->input[0] );

	vtkImageMathematics * sum = vtkImageMathematics::New();
	sum->SetOperationToAdd();
	sum->SetInput1( average );

	// build initial estimate
	this->metric->setInput1( average );
	
	vtkImageMathematics * multiplyByK = vtkImageMathematics::New();
	multiplyByK->SetOperationToMultiplyByK();
	multiplyByK->SetInput1( sum->GetOutput() );
	multiplyByK->SetConstantK( .5 );

	for ( int k = 1; k < this->input.size(); k++ ){
		// compute alignment
		this->metric->setInput2( this->input[k] );
		this->metric->getDistance();
		// overwrite input with rotated version
		this->input[k]->DeepCopy( this->metric->getAligned() );
		sum->SetInput2( this->input[k] );
		sum->Update();
		multiplyByK->Update();
		average->ShallowCopy( multiplyByK->GetOutput() );
	}

	vtkImageMathematics * substract = vtkImageMathematics::New();
	substract->SetOperationToSubtract();
	substract->SetInput1( average );

	vtkImageMathematics * multiplyByK2 = vtkImageMathematics::New();
	multiplyByK2->SetOperationToMultiplyByK();
	multiplyByK2->SetConstantK( 1.0 / ( this->input.size() + 0.0 ) );

	// refinement loop
	for ( int i = 0; i < numberOfIterations; i++ ){
		for ( int j = 1; j < this->input.size(); j++ ){
			// substract input[i] / size(input)
			multiplyByK2->SetInput1( this->input[j] );
			multiplyByK2->Update();
			substract->SetInput2( multiplyByK2->GetOutput() );
			substract->Update();
			this->metric->setInput1( substract->GetOutput() );
			this->metric->setInput2( this->input[j] );
			this->metric->getDistance();
			sum->SetInput1( this->metric->getAligned() );
			sum->SetInput2( substract->GetOutput() );
			sum->Update();
			multiplyByK->SetInput1( sum->GetOutput() );
			multiplyByK->Update();
			average->ShallowCopy( multiplyByK->GetOutput() );
		}
	}
	output->DeepCopy( average );
	average->Delete();
}