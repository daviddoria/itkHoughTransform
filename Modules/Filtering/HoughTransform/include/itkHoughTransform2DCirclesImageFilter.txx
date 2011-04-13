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
#ifndef __itkHoughTransform2DCirclesImageFilter_txx
#define __itkHoughTransform2DCirclesImageFilter_txx

#include "itkHoughTransform2DCirclesImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"

namespace itk
{
template< typename TInputImageType >
HoughTransform2DCirclesImageFilter< TInputImageType >
::HoughTransform2DCirclesImageFilter()
{
  m_Threshold = 0;
  m_MinimumRadius = 0;  // by default
  m_MaximumRadius = 10; // by default
  m_SigmaGradient = 1;  // Scale of the DoG filter
  m_DiscRadiusRatio = 1;
  m_Variance = 10;
  m_OldModifiedTime = 0;
  m_OldNumberOfCircles = 0;
  m_SweepAngle = 0.0;
  m_NumberOfCircles = 1;
}

template< typename TInputImageType >
void
HoughTransform2DCirclesImageFilter< TInputImageType >
::SetRadius(double radius)
{
  this->SetMinimumRadius(radius);
  this->SetMaximumRadius(radius);
}


/** Get the list of circles. This recomputes the circles */
template< typename TInputPixelType, typename TOutputPixelType >
typename HoughTransform2DCirclesImageFilter< TInputPixelType, TOutputPixelType >::CirclesListType &
HoughTransform2DCirclesImageFilter< TInputPixelType, TOutputPixelType >
::GetCircles(unsigned int n)
{
  if ( ( this->GetMTime() == m_OldModifiedTime ) && ( n == m_OldNumberOfCircles ) )
    {
    // if the filter has not been updated
    return m_CirclesList;
    }

  m_CirclesList.clear();

  BlurAccumulator();

  Index< 2 > index;

  CirclesListSizeType circles = 0;
  bool         found;

  // Find maxima
  do
    {
    minMaxCalculator->SetImage(postProcessImage);
    minMaxCalculator->ComputeMaximum();
    InternalImageType::PixelType max = minMaxCalculator->GetMaximum();

    found = false;
    for ( it_input.GoToBegin(); !it_input.IsAtEnd(); ++it_input )
      {
      if ( it_input.Get() == max )
        {
        // Create a Line Spatial Object
        CirclePointer Circle = CircleType::New();
        Circle->SetId(circles);
        Circle->SetRadius( m_RadiusImage->GetPixel( it_input.GetIndex() ) );

        CircleType::VectorType center;
        center[0] = it_input.GetIndex()[0];
        center[1] = it_input.GetIndex()[1];
        Circle->GetObjectToParentTransform()->SetOffset(center);
        Circle->ComputeBoundingBox();

        m_CirclesList.push_back(Circle);

        // Remove a black disc from the hough space domain
        for ( double angle = 0; angle <= 2 * itk::Math::pi; angle += itk::Math::pi / 1000 )
          {
          for ( double length = 0; length < m_DiscRadiusRatio * Circle->GetRadius()[0]; length += 1 )
            {
            index[0] = (IndexValueType)( it_input.GetIndex()[0] + length * vcl_cos(angle) );
            index[1] = (IndexValueType)( it_input.GetIndex()[1] + length * vcl_sin(angle) );
            if ( postProcessImage->GetLargestPossibleRegion().IsInside(index) )
              {
              postProcessImage->SetPixel(index, 0);
              }
            }
          }
        minMaxCalculator->SetImage(postProcessImage);
        minMaxCalculator->ComputeMaximum();
        max = minMaxCalculator->GetMaximum();

        circles++;
        found = true;
        if ( circles == m_NumberOfCircles ) { break; }
        }
      }
    }
  while ( ( circles < m_NumberOfCircles ) && ( found ) );

  m_OldModifiedTime = this->GetMTime();
  m_OldNumberOfCircles = m_CirclesList.size();
  return m_CirclesList;
}

/** Print Self information */
template< typename TInputImageType >
void
HoughTransform2DCirclesImageFilter< TInputImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  //os << "Threshold: " << m_Threshold << std::endl;

}
} // end namespace

#endif
