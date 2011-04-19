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
template<typename TModelParameter>
HoughTransform2DLinesImageFilter<TModelParameter>
::HoughTransform2DLinesImageFilter()
{

}

/** All model paramters except one are held constant and an input point is given. The remaining model paramter is solved. */
template<typename TModelParameter>
float HoughTransform2DLinesImageFilter<TModelParameter>
::SolveModel(itk::Point<2> point, itk::FixedArray<float, VModelDimension> parameters, unsigned int parameterToSolve)
{
  //Dimension 0 is R
  //Dimension 1 is Theta
  
  if(parameterToSolve == 1)
    {
    std::cerr << "Cannot solve this paramter." << std::endl;
    return 0;
    }
  
  if(parameterToSolve == 0)
    {
    //R = x*vcl_cos(Theta)+y*vcl_sin(Theta);
    return point[0] * vcl_cos(parameters[1]) + points[1] * vcl_sin(parameters[1]);
    }
}

template<typename TModelParameter>
LineSpatialObject::Pointer HoughTransform2DLinesImageFilter<TModelParameter>
::CreateObject(itk::FixedArray<float, VModelDimension> parameters)
{
  // A LineSpatialObject is just a list of points (2 in this case). We intersect the line with the bounding box of the data
  // to choose the two points to create the line with.

/*
  // Line equation is y = (- cos(Theta)/sin(Theta))x + (r/sin(Theta))
  // We can find two points on this line by computing y at x=0 and x=1. This will give us a ray.
  
  //float y0 = -vcl_cos(Theta)/vcl_sin(Theta) * 0.0 + R/vcl_sin(Theta);
  float y0 = parameters[0]/vcl_sin(parameters[1]);
  
  //float y1 = -vcl_cos(Theta)/vcl_sin(Theta) * 0.0 + R/vcl_sin(Theta);
  float y1 = -vcl_cos(parameters[1])/vcl_sin(parameters[1]) * 1.0 + parameters[0]/vcl_sin(parameters[1]);
*/

  // Find bounding box
  pointSet->ComputeBoundingBox();
  const PointSetType::BoundingBoxType* boundingBox = pointSet->GetBoundingBox();

  // Find intersections of ray with bounding box parallel to 0 axis
  LineType::LinePointType intersection1;
  
  // First side of bounding box
  float y0 = -vcl_cos(Theta)/vcl_sin(Theta) * boundingBox->GetBounds()[0] + R/vcl_sin(Theta);
  if(y0 > boundingBox->GetBounds()[1] && y0 < boundingBox->GetBounds()[2])
    {
    intersection1.SetPosition(boundingBox->GetBounds()[0], y0);
    }
  
  // Second side of bounding box

...

  LineType::LinePointType intersection2;
  point.SetPosition(10,i);

  std::vector<LineType::LinePointType> points(2);
  points[0] = intersection1;
  points[1] = intersection2;

  // Create the line
  LineType::Pointer line = LineType::New();
  line->SetPoints(points);
  
  return line;
}

/** Print Self information */
template<typename TModelParameter>
void
HoughTransform2DLinesImageFilter<TModelParameter>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

}
} // end namespace

#endif
