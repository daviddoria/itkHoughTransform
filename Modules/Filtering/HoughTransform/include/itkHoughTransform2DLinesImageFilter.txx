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
HoughTransform2DLinesImageFilter
::HoughTransform2DLinesImageFilter()
{

}

/** All model paramters except one are held constant and an input point is given. The remaining model paramter is solved. */
float HoughTransform2DLinesImageFilter
::SolveModel(itk::FixedArray<float, VModelDimension> parameters, unsigned int parameterToSolve)
{
  // R = x*vcl_cos(Theta)+y*vcl_sin(Theta)
  
  // x = (R - y*vcl_sin(Theta))/vcl_cos(Theta)
  
  // y = (R - x*vcl_cos(Theta))/vcl_sin(Theta)

}

LineSpatialObject::Pointer HoughTransform2DLinesImageFilter::CreateObject(itk::FixedArray<float, VModelDimension> parameters)
{
  // A LineSpatialObject is just a list of points (2 in this case). How do we choose the two points to create the line with if we
  // only know the perpendicular distance from the origin (R) and the angle (Theta)?
}

/** Print Self information */
template< typename TInputPixelType, typename TOutputPixelType >
void
HoughTransform2DLinesImageFilter< TInputPixelType, TOutputPixelType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

}
} // end namespace

#endif
