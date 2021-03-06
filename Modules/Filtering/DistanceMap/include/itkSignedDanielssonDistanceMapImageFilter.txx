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
#ifndef __itkSignedDanielssonDistanceMapImageFilter_txx
#define __itkSignedDanielssonDistanceMapImageFilter_txx

#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
/**
 *    Constructor
 */
template< class TInputImage, class TOutputImage >
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::SignedDanielssonDistanceMapImageFilter()
{
  this->SetNumberOfRequiredOutputs(3);

  OutputImagePointer voronoiMap = OutputImageType::New();
  this->SetNthOutput( 1, voronoiMap.GetPointer() );

  VectorImagePointer distanceVectors = VectorImageType::New();
  this->SetNthOutput( 2, distanceVectors.GetPointer() );

  //Default values
  this->m_SquaredDistance     = false;  //Should we remove this ?
                                        //doesn't make sense in a SignedDaniel
  this->m_UseImageSpacing     = false;
  this->m_InsideIsPositive    = false;
}

/** This is overloaded to create the VectorDistanceMap output image */
template< class TInputImage, class TOutputImage >
typename SignedDanielssonDistanceMapImageFilter<
  TInputImage, TOutputImage >::DataObjectPointer
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::MakeOutput(unsigned int idx)
{
  if  ( idx == 2 )
    {
    return static_cast< DataObject * >( VectorImageType::New().GetPointer() );
    }
  return Superclass::MakeOutput(idx);
}

/**
 *  Return the distance map Image pointer
 */
template< class TInputImage, class TOutputImage >
typename SignedDanielssonDistanceMapImageFilter<
  TInputImage, TOutputImage >::OutputImageType *
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::GetDistanceMap(void)
{
  return dynamic_cast< OutputImageType * >(
           this->ProcessObject::GetOutput(0) );
}

/**
 *  Return Closest Points Map
 */
template< class TInputImage, class TOutputImage >
typename
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >::OutputImageType *
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::GetVoronoiMap(void)
{
  return dynamic_cast< OutputImageType * >(
           this->ProcessObject::GetOutput(1) );
}

/**
 *  Return the distance vectors
 */
template< class TInputImage, class TOutputImage >
typename SignedDanielssonDistanceMapImageFilter<
  TInputImage, TOutputImage >::VectorImageType *
SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::GetVectorDistanceMap(void)
{
  return dynamic_cast< VectorImageType * >(
           this->ProcessObject::GetOutput(2) );
}

/**
 *  Compute Distance and Voronoi maps by calling
 * DanielssonDistanceMapImageFilter twice.
 */
template< class TInputImage, class TOutputImage >
void SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  //Set up mini pipeline filter
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef DanielssonDistanceMapImageFilter<
    InputImageType, OutputImageType >  FilterType;
  typename FilterType::Pointer filter1 = FilterType::New();
  typename FilterType::Pointer filter2 = FilterType::New();

  filter1->SetInputIsBinary(true);    // Force signed distance map to work on
  filter2->SetInputIsBinary(true);    // binary images
  filter1->SetUseImageSpacing(m_UseImageSpacing);
  filter2->SetUseImageSpacing(m_UseImageSpacing);
  filter1->SetSquaredDistance(m_SquaredDistance);
  filter2->SetSquaredDistance(m_SquaredDistance);

  //Invert input image for second Danielsson filter
  typedef typename InputImageType::PixelType                InputPixelType;
  typedef Functor::InvertIntensityFunctor< InputPixelType > FunctorType;
  typedef UnaryFunctorImageFilter< InputImageType,
                                   InputImageType,
                                   FunctorType >             InverterType;

  typename InverterType::Pointer inverter = InverterType::New();

  inverter->SetInput( this->GetInput() );

  //Dilate the inverted image by 1 pixel to give it the same boundary
  //as the univerted input.

  typedef BinaryBallStructuringElement<
    InputPixelType,
    InputImageDimension  > StructuringElementType;

  typedef BinaryDilateImageFilter<
    InputImageType,
    InputImageType,
    StructuringElementType >     DilatorType;

  typename DilatorType::Pointer dilator = DilatorType::New();

  StructuringElementType structuringElement;
  structuringElement.SetRadius(1);    // 3x3 structuring element
  structuringElement.CreateStructuringElement();
  dilator->SetKernel(structuringElement);
  dilator->SetDilateValue(1);

  filter1->SetInput( this->GetInput() );
  dilator->SetInput( inverter->GetOutput() );
  filter2->SetInput( dilator->GetOutput() );

  //Subtract Distance maps results of the two Danielsson filters
  typedef SubtractImageFilter< OutputImageType, OutputImageType,
                               OutputImageType > SubtracterType;

  typename SubtracterType::Pointer subtracter = SubtracterType::New();

  if ( m_InsideIsPositive )
    {
    subtracter->SetInput1( filter2->GetDistanceMap() );
    subtracter->SetInput2( filter1->GetDistanceMap() );
    }
  else
    {
    subtracter->SetInput2( filter2->GetDistanceMap() );
    subtracter->SetInput1( filter1->GetDistanceMap() );
    }

  subtracter->Update();
  filter1->Update();
  filter2->Update();

  // Register progress
  progress->RegisterInternalFilter(filter1, .5f);

  // Graft outputs
  this->GraftNthOutput( 0, subtracter->GetOutput() );

  // we must use not use this->GetOutput method because the
  // output's are of different types then ImageSource
  this->GraftNthOutput( 1, filter1->GetVoronoiMap() );
  this->GraftNthOutput( 2, filter1->GetVectorDistanceMap() );
}

// end GenerateData()

/**
 *  Print Self
 */
template< class TInputImage, class TOutputImage >
void SignedDanielssonDistanceMapImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Signed Danielson Distance: " << std::endl;
  os << indent << "Use Image Spacing : " << m_UseImageSpacing << std::endl;
  os << indent << "Squared Distance  : " << m_SquaredDistance << std::endl;
  os << indent << "Inside is positive  : " << m_InsideIsPositive << std::endl;
}
} // end namespace itk

#endif
