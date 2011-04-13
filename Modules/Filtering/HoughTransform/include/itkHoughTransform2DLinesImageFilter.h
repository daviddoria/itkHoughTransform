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
 * This filter derives from ImageToImageFilter
 * The input is an image, and all pixels above some threshold are those
 * to be extracted. The output is the image of the accumulator.
 * GetLines() returns a list of LinesSpatialObjects
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
 * The size of the array depends on the AngleAxisSize that could be set
 * (500 by default) for the angle axis. The distance axis depends on the
 * size of the diagonal of the input image.
 *
 * \ingroup ImageFeatureExtraction
 * \sa LineSpatialObject
 *
 * */

template< typename TInputImageType>
class ITK_EXPORT HoughTransform2DLinesImageFilter:
  public HoughTransform<2>
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
  /*
  typedef typename LineType::Pointer LinePointer;
  typedef std::list< LinePointer >   LinesListType;
  typedef LineType::LinePointType    LinePointType;
  */

  typedef typename LinesListType::size_type LinesListSizeType;

  /** Standard "Superclass" typedef. */
  typedef HoughTransform<2> Superclass;

  /** Image index typedef */
  typedef typename InputImageType::IndexType IndexType;

  /** Image pixel value typedef */
  typedef typename InputImageType::PixelType PixelType;

  /** Typedef to describe the output image region type. */
  typedef typename InputImageType::RegionType OutputImageRegionType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(HoughTransform2DLinesImageFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Method for evaluating the implicit function over the image. */
  void GenerateData();

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToOutputCheck,
                   ( Concept::Convertible< int, TOutputPixelType > ) );
  itkConceptMacro( InputGreaterThanFloatCheck,
                   ( Concept::GreaterThanComparable< PixelType, float > ) );
  itkConceptMacro( OutputPlusIntCheck,
                   ( Concept::AdditiveOperators< TOutputPixelType, int > ) );
  /** End concept checking */
#endif

protected:

  HoughTransform2DLinesImageFilter();
  virtual ~HoughTransform2DLinesImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const;

private:

  HoughTransform2DLinesImageFilter(const Self &);
  void operator=(const Self &);

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHoughTransform2DLinesImageFilter.txx"
#endif

#endif
