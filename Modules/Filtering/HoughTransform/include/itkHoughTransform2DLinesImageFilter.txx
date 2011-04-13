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
#ifndef __itkHoughTransform2DLinesImageFilter_txx
#define __itkHoughTransform2DLinesImageFilter_txx

#include "itkHoughTransform2DLinesImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCastImageFilter.h"

namespace itk
{
/** Constructor */
template< typename TInputImageType >
HoughTransform2DLinesImageFilter< TInputPixelType >
::HoughTransform2DLinesImageFilter()
{

}

/** All model paramters except one are held constant and an input point is given. The remaining model paramter is solved. */
template< unsigned int VModelDimension >
float
HoughTransform< VModelDimension >
::SolveModel(itk::FixedArray<float, VModelDimension> parameters, unsigned int parameterToSolve)
{
  // R = x*vcl_cos(Theta)+y*vcl_sin(Theta)

}


/** Get the list of lines. This recomputes the lines */
template< typename TInputPixelType, typename TOutputPixelType >
typename HoughTransform2DLinesImageFilter< TInputPixelType, TOutputPixelType >::LinesListType &
HoughTransform2DLinesImageFilter< TInputPixelType, TOutputPixelType >
::GetLines(unsigned int n)
{
  // if the filter has not been updated
  if ( ( this->GetMTime() == m_OldModifiedTime ) && ( n == m_OldNumberOfLines ) )
    {
    return m_LinesList;
    }

  m_LinesList.clear();

  /** Blur the accumulator in order to find the maximum */
  typedef float                              InternalImagePixelType;
  typedef Image< InternalImagePixelType, 2 > InternalImageType;

  OutputImagePointer outputImage = this->GetOutput(0);

  if ( !outputImage )
    {
    itkExceptionMacro("Update() must be called before GetLines().");
    }

  /** xxxConvert the accumulator output image type to internal image type */
  typedef CastImageFilter< OutputImageType, InternalImageType > CastImageFilterType;

  typename CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
  castImageFilter->SetInput(outputImage);

  typedef DiscreteGaussianImageFilter< InternalImageType, InternalImageType > GaussianFilterType;
  typename GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();

  // the output is the accumulator image
  gaussianFilter->SetInput( castImageFilter->GetOutput() );
  double variance[2];
  variance[0] = m_Variance;
  variance[1] = m_Variance;
  gaussianFilter->SetVariance(variance);
  gaussianFilter->Update();
  InternalImageType::Pointer postProcessImage = gaussianFilter->GetOutput();

  typedef MinimumMaximumImageCalculator< InternalImageType > MinMaxCalculatorType;
  typename MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();
  itk::ImageRegionIterator< InternalImageType >
  it_input( postProcessImage, postProcessImage->GetLargestPossibleRegion() );

  itk::Index< 2 > index;

  unsigned int lines = 0;
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
        // Create the line
        LineType::PointListType list; // insert two points per line

        double radius = it_input.GetIndex()[0];
        double teta   = ( ( it_input.GetIndex()[1] ) * 2 * itk::Math::pi / this->GetAngleResolution() ) - itk::Math::pi;
        double Vx = radius * vcl_cos(teta);
        double Vy = radius * vcl_sin(teta);
        double norm = vcl_sqrt(Vx * Vx + Vy * Vy);
        double VxNorm = Vx / norm;
        double VyNorm = Vy / norm;

        if ( ( teta <= 0 ) || ( teta >= itk::Math::pi / 2 ) )
          {
          if ( teta >= itk::Math::pi / 2 )
            {
            VyNorm = -VyNorm;
            VxNorm = -VxNorm;
            }

          LinePointType p;
          p.SetPosition(Vx, Vy);
          list.push_back(p);
          p.SetPosition(Vx - VyNorm * 5, Vy + VxNorm * 5);
          list.push_back(p);
          }
        else // if teta>0
          {
          LinePointType p;
          p.SetPosition(Vx, Vy);
          list.push_back(p);
          p.SetPosition(Vx - VyNorm * 5, Vy + VxNorm * 5);
          list.push_back(p);
          } // end if(teta>0)

        // Create a Line Spatial Object
        LinePointer Line = LineType::New();
        Line->SetId(lines);
        Line->SetPoints(list);
        Line->ComputeBoundingBox();

        m_LinesList.push_back(Line);

        // Remove a black disc from the hough space domain
        for ( double angle = 0; angle <= 2 * itk::Math::pi; angle += itk::Math::pi / 1000 )
          {
          for ( double length = 0; length < m_DiscRadius; length += 1 )
            {
            index[0] = (IndexValueType)( it_input.GetIndex()[0] + length * vcl_cos(angle) );
            index[1] = (IndexValueType)( it_input.GetIndex()[1] + length * vcl_sin(angle) );
            if ( postProcessImage->GetBufferedRegion().IsInside(index) )
              {
              postProcessImage->SetPixel(index, 0);
              }
            }
          }
        minMaxCalculator->SetImage(postProcessImage);
        minMaxCalculator->ComputeMaximum();
        max = minMaxCalculator->GetMaximum();

        lines++;
        found = true;
        if ( lines == m_NumberOfLines ) { break; }
        }
      }
    }
  while ( ( lines < m_NumberOfLines ) && ( found ) );

  m_OldModifiedTime = this->GetMTime();
  m_OldNumberOfLines = m_LinesList.size();
  return m_LinesList;
}

/** Print Self information */
template< typename TInputPixelType, typename TOutputPixelType >
void
HoughTransform2DLinesImageFilter< TInputPixelType, TOutputPixelType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  //os << "Threshold: " << m_Threshold << std::endl;

}
} // end namespace

#endif
