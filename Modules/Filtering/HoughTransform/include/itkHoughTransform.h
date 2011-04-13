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

#include "itkImageSource.h"

namespace itk
{
/**
 * \class HoughTransform
 * \brief An abstract base class which sets up the structure for a Hough
 * Transform to detect a particular object.
 *
 * This class provides the methods which produce the Hough accumulator array.
 *
 * Subclasses may either operate on an itk::Image or an itk::PointSet.
 *
 * The output of this class is the Hough accumulator array.
 * \ingroup ImageFeatureExtraction
 *
 * */

template<unsigned int VModelDimension >
class ITK_EXPORT HoughTransform:
  public ImageSource< itk::Image< float, VModelDimension > >
{
public:

  /** Standard "Self" typedef. */
  typedef HoughTransform Self;

  /** Output Image typedef */
  typedef Image< float, VModelDimension >      OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  /** Smart pointer typedef support. */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Standard "Superclass" typedef. */
  typedef ImageSource< OutputImageType > Superclass;

  /** Image index typedef */
  typedef typename OutputImageType::IndexType IndexType;

  /** Image pixel value typedef */
  typedef typename OutputImageType::PixelType PixelType;

  /** Typedef to describe the output image region type. */
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** Method for evaluating the implicit function over the image. */
  void GenerateData();

  /** Set the threshold above which the filter should consider
      the point as a valid point */
  itkSetMacro(Threshold, float);

  /** Get the threshold value */
  itkGetConstMacro(Threshold, float);

  /** Set the number of bins in each dimension of the accumulator array */
  itkSetMacro(AccumulatorArrayDimensions, itk::FixedArray<unsigned int, VModelDimension>);

  /** Get the number of bins in each dimension of the accumulator array */
  itkGetConstMacro(AccumulatorArrayDimensions, itk::FixedArray<unsigned int, VModelDimension>);

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

  /** Blur the accumulator array. */
  void BlurAccumulator();

  /** All model paramters except one are held constant and an input point is given. The remaining model paramter is solved. */
  virtual float SolveModel(itk::FixedArray<float, VModelDimension> parameters, unsigned int parameterToSolve) = 0;

  /** Object typedef */
  typedef SpatialObject< 2 >    	ObjectType;
  typedef typename ObjectType::Pointer 	ObjectPointer;
  typedef std::list< ObjectPointer >   	ObjectListType;

  typedef typename ObjectListType::size_type ObjectListSizeType;

  virtual ObjectListType & GetObjects() = 0;

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

  /** The array of the number of bins in each dimension of the Hough accumulator array */
  itk::FixedArray<unsigned int, VModelDimension> m_AccumulatorArrayDimensions;

  /** The array of the min/max pairs of each dimension of the Hough accumulator array */
  itk::FixedArray<std::pair<float, float>, VModelDimension> m_AccumulatorArrayBounds;

  /** A vector of all valid combinations of the paramters */
  std::vector<itk::FixedArray<float, VModelDimension> > m_ParameterList;

  /** The threshold above which the filter should consider
      the point as a valid point*/
  float              m_Threshold;

  /** The simplified accumulator */
  OutputImagePointer m_SimplifyAccumulator;

  /** The variance of the Gaussian kernel used to smooth the accumulator */
  float              m_Variance;

};
} // end namespace itk

//#ifndef ITK_MANUAL_INSTANTIATION
//#include "itkHoughTransform.txx"
//#endif

#endif
