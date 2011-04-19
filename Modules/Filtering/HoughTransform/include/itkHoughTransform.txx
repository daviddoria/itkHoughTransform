/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkHoughTransform_txx
#define __itkHoughTransform_txx

#include "itkHoughTransform.h"

namespace itk
{
/** Constructor */
template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::HoughTransform()
{
  m_NumberOfObjectsToFind = 1;
  m_Variance = 5;
}

template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  // call the superclass' implementation of this method
  Superclass::EnlargeOutputRequestedRegion(output);

  output->SetRequestedRegionToLargestPossibleRegion();
}

template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer     output = this->GetOutput();

  if ( !input || !output )
    {
    return;
    }

  // Compute the size of the output image
  typename InputImageType::RegionType region;
  Size< VModelDimension > size;

  for(unsigned int dimension = 0; dimension < VModelDimension; ++dimension)
    {
    size[dimension] = m_AccumulatorArrayDimensions[dimension];
    }

  region.SetSize(size);
  region.SetIndex( input->GetLargestPossibleRegion().GetIndex() );

  output->SetLargestPossibleRegion(region);
}


/** Blur the accumulator image */
template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::BlurAccumulator()
{
  typedef Image< float, VModelDimension > InternalImageType;

  OutputImagePointer outputImage = OutputImageType::New();
  outputImage->SetRegions( this->GetOutput(0)->GetLargestPossibleRegion() );
  outputImage->SetOrigin( this->GetOutput(0)->GetOrigin() );
  outputImage->SetSpacing( this->GetOutput(0)->GetSpacing() );
  outputImage->SetDirection( this->GetOutput(0)->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer(0);

  ImageRegionConstIteratorWithIndex< OutputImageType > image_it( this->GetOutput(0),  this->GetOutput(
                                                                   0)->GetRequestedRegion() );
  image_it.GoToBegin();

  ImageRegionIterator< InternalImageType > it( outputImage,  outputImage->GetRequestedRegion() );

  while ( !image_it.IsAtEnd() )
    {
    it.Set( image_it.Get() );
    ++image_it;
    ++it;
    }

  typedef DiscreteGaussianImageFilter< OutputImageType, InternalImageType > GaussianFilterType;
  typename GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();

  gaussianFilter->SetInput(outputImage); // the output is the accumulator image
  double variance[2];
  variance[0] = m_Variance;
  variance[1] = m_Variance;
  gaussianFilter->SetVariance(variance);
  gaussianFilter->Update();

}

/** Generate and blur the accumulator image. */
template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::GenerateData()
{
  // Get the output pointer
  OutputImagePointer     outputImage = this->GetOutput(0);

  // Allocate the output
  outputImage->FillBuffer(0);

  // For every point, solve all of the equations corresponding to every voxel in the accumulator array
  itk::ImageRegionIterator<ImageType> iterator(outputImage, outputImage->GetLargestPossibleRegion());
  
  for(unsigned int pointId = 0; pointId < numberOfPoints; ++pointId)
    {
    while ( !iterator.IsAtEnd() )
      {
      itk::Index<VModelDimension> accumulatorIndex = iterator.GetIndex();
    
      // The real work is dispatched to a subclass
      SolveModel(accumulatorIndex);

      outputImage->SetPixel(index, outputImage->GetPixel(index) + 1);
      ++iterator;
      } // end while over accumulator array
    } // end for over points
    
  BlurAccumulator();
}

/** Get the list of NumberOfObjects objects from the accumulator array. */
template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
typename HoughTransform< VModelDimension>::ObjectListType &
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::GetObjects(unsigned int n)
{
  ObjectListType objects;

  typedef MinimumMaximumImageCalculator< InternalImageType > MinMaxCalculatorType;
  typename MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();

  // Find maxima
  unsigned int numberOfMaximaFound = 0;
  while(numberOfMaximaFound < this->m_NumberOfObjectsToFind)
    {
    minMaxCalculator->SetImage(postProcessImage);
    minMaxCalculator->ComputeMaximum();

    itk::Index<VModelDimension> indexOfMaximum = minMaxCalculator->GetIndexOfMaxmum();

    // Create the object
    ObjectPointer object = ObjectType::New();
    CreateObject(indexOfMaximum, object);
  
    objects.push_back(object);

    ClearRegion(indexOfMaximum);
    
    numberOfMaximaFound++;
    }

  return objects;
}

template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::ClearRegion(itk::Index<VModelDimension> index)
{
  // Zero a sphere around the specified index
  
  typedef itk::BinaryBallStructuringElement<2> StructuringElementType;
  StructuringElementType::RadiusType elementRadius;
  elementRadius.Fill(m_Variance);
  StructuringElementType structuringElement =
    StructuringElementType::Box(elementRadius);

  typedef itk::ShapedNeighborhoodIterator<ImageType> IteratorType;
  IteratorType iterator(radius, image, image->GetLargestPossibleRegion());
  iterator.SetNeighborhood(structuringElement);
}

/** Print Self information */
template<typename TPoint, unsigned int VInputDimension, typename TModelParameter, unsigned int VModelDimension, typename TSpatialObject >
void
HoughTransform<TPoint, VInputDimension, TModelParameter, VModelDimension, TSpatialObject >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Number of objects to find: " << m_NumberOfObjectsToFind << std::endl;
  os << "Accumulator blur variance: " << m_Variance << std::endl;

}
} // end namespace

#endif
