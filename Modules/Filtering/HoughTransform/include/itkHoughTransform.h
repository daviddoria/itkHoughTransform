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
#ifndef __itkHoughTransform_h
#define __itkHoughTransform_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkImageToImageFilter.h"
#include "itkLineSpatialObject.h"

namespace itk
{
/**
 * \class HoughTransform
 * \brief An abstract base class which sets up the structure for a Hough
 * Transform to detect a particular object.
 *
 * \ingroup ImageFeatureExtraction
 *
 * */

template< typename TInputPixelType, typename TOutputPixelType, unsigned int VModelDimension >
class ITK_EXPORT HoughTransform:
  public ImageToImageFilter< Image< TInputPixelType, 2 >, Image< TOutputPixelType, 2 > >
{
public:

  /** Standard "Self" typedef. */
  typedef HoughTransform Self;

  /** Input Image typedef */
  typedef Image< TInputPixelType, 2 >           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  /** Output Image typedef */
  typedef Image< TOutputPixelType, 2 >      OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  /** Smart pointer typedef support. */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Standard "Superclass" typedef. */
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;

  /** Image index typedef */
  typedef typename InputImageType::IndexType IndexType;

  /** Image pixel value typedef */
  typedef typename InputImageType::PixelType PixelType;

  /** Typedef to describe the output image region type. */
  typedef typename InputImageType::RegionType OutputImageRegionType;

  /** Method for evaluating the implicit function over the image. */
  void GenerateData();

  /** Set the threshold above which the filter should consider
      the point as a valid point */
  itkSetMacro(Threshold, float);

  /** Get the threshold value */
  itkGetConstMacro(Threshold, float);

  /** Set the resolution of Hough dimension 0 */
  itkSetMacro(Dimension0Resolution, float);

  /** Get the resolution of Hough dimension 0 */
  itkGetConstMacro(Dimension0Resolution, float);

  /** Set the resolution of Hough dimension 1 */
  itkSetMacro(Dimension1Resolution, float);

  /** Get the resolution of Hough dimension 1 */
  itkGetConstMacro(Dimension1Resolution, float);

  /** Set the resolution of Hough dimension 2 */
  itkSetMacro(Dimension2Resolution, float);

  /** Get the resolution of Hough dimension 2 */
  itkGetConstMacro(Dimension2Resolution, float);

  /** Simplify the accumulator */
  void Simplify(void);

  /** Get the Simplified accumulator */
  itkGetObjectMacro(SimplifyAccumulator, OutputImageType);

  /** Set/Get the number of objects to extract */
  itkSetMacro(NumberOfObjects, unsigned int);
  itkGetConstMacro(NumberOfObjects, unsigned int);

  /** Set the variance of the gaussian bluring for the accumulator */
  itkSetMacro(Variance, float);
  itkGetConstMacro(Variance, float);

  /** Object typedef */
  typedef SpatialObject< 2 >    	ObjectType;
  typedef typename ObjectType::Pointer 	ObjectPointer;
  typedef std::list< ObjectPointer >   	ObjectListType;

  typedef typename ObjectListType::size_type ObjectListSizeType;

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

  HoughTransform();
  virtual ~HoughTransform() {}

  void PrintSelf(std::ostream & os, Indent indent) const;

  /** HoughTransform needs the entire input. Therefore
   * it must provide an implementation GenerateInputRequestedRegion().
   * \sa ProcessObject::GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** HoughTransform's output is the accumulator
   * array.  The size of the output is a function of the size of the
   * input. Since this output is a different
   * size than the input, it must provide an implementation of
   * GenerateOutputInformation.
   * \sa ProcessObject::GenerateOutputRequestedRegion() */
  void GenerateOutputInformation();

  /** HoughTransform must produce the entire output */
  void EnlargeOutputRequestedRegion(DataObject *output);

private:

  HoughTransform(const Self &);
  void operator=(const Self &);

  float              m_AxesResolution[VModelDimension];
  float              m_Dimension1Resolution;
  float              m_Dimension2Resolution;
  float              m_Threshold;
  OutputImagePointer m_SimplifyAccumulator;
  float              m_Variance;
  unsigned long      m_OldModifiedTime;

};
} // end namespace itk

//#ifndef ITK_MANUAL_INSTANTIATION
//#include "itkHoughTransform.txx"
//#endif

#endif
