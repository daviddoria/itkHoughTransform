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
#ifndef __itkHoughTransform2DLinesImageFilter_h
#define __itkHoughTransform2DLinesImageFilter_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkImageToImageFilter.h"
#include "itkLineSpatialObject.h"

namespace itk
{
/**
 * \class HoughTransform2DLinesImageFilter
 * \brief Performs the Hough Transform to find 2D straight lines
 *        in a 2D image.
 *
 * This filter derives from HoughTransform
 * The input is a PointSet. If you have an image, you must threshold
 * and extract points.
 * 
 * The output is the image of the accumulator array.
 * GetObjects() returns a list of LineSpatialObjects.
 *
 * Lines are parameterized in the form: R = x*vcl_cos(Theta)+y*vcl_sin(Theta)
 * where R is the perpendicular distance from the origin and Theta
 * the angle with the normal.
 *
 * The output is the accumulator array:
 *    -Dimension 0 represents the distance R from the corner
 *     to the line
 *    -Dimension 1 represents the angle between the X axis
 *     and the normal to the line.
 *
 * \ingroup ImageFeatureExtraction
 * \sa LineSpatialObject
 *
 *
 * \ingroup ITK-ImageFeature
 * \wikiexample{Conversions/HoughTransform2DLinesImageFilter,HoughTransform2DLinesImageFilter}
 */

class ITK_EXPORT HoughTransform2DLinesImageFilter:
  public HoughTransform<float, 2, LineSpatialObject>
{
public:

  /** Standard "Self" typedef. */
  typedef HoughTransform2DLinesImageFilter Self;

  /** Input Image typedef */
  typedef Image< TInputPixelType, 2 >           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  /** Smart pointer typedef support. */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Line typedef */
  typedef LineSpatialObject< 2 >     ObjectType;
  typedef typename LineType::Pointer LinePointer;

  /** Standard "Superclass" typedef. */
  typedef HoughTransform<2> Superclass;

  /** Run-time type information (and related methods). */
  itkTypeMacro(HoughTransform2DLinesImageFilter, HoughTransform);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Compute the entire accumulator array. */
  void GenerateData();

protected:

  HoughTransform2DLinesImageFilter();
  virtual ~HoughTransform2DLinesImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const;

  virtual ObjectPointer CreateObject(itk::FixedArray<float, VModelDimension>);
  
private:

  HoughTransform2DLinesImageFilter(const Self &);
  void operator=(const Self &);

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHoughTransform2DLinesImageFilter.txx"
#endif

#endif
