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
#ifndef __itkShapeRelabelLabelMapFilter_h
#define __itkShapeRelabelLabelMapFilter_h

#include "itkInPlaceLabelMapFilter.h"
#include "itkLabelObject.h"
#include "itkShapeLabelObjectAccessors.h"

namespace itk
{
/** \class ShapeRelabelLabelMapFilter
 * \brief Relabels objects according to their shape attributes.
 *
 * The ShapeRelabelImageFilter relabels a label collection image according to the shape attributes of
 * the objects. The label produced are always consecutives.
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/584  or
 * http://www.insight-journal.org/browse/publication/176
 *
 * \author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \sa ShapeLabelObject, RelabelComponentImageFilter
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 * \ingroup ITK-Review
 */
template< class TImage >
class ITK_EXPORT ShapeRelabelLabelMapFilter:
  public InPlaceLabelMapFilter< TImage >
{
public:
  /** Standard class typedefs. */
  typedef ShapeRelabelLabelMapFilter      Self;
  typedef InPlaceLabelMapFilter< TImage > Superclass;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;

  /** Some convenient typedefs. */
  typedef TImage                              ImageType;
  typedef typename ImageType::Pointer         ImagePointer;
  typedef typename ImageType::ConstPointer    ImageConstPointer;
  typedef typename ImageType::PixelType       PixelType;
  typedef typename ImageType::IndexType       IndexType;
  typedef typename ImageType::LabelObjectType LabelObjectType;

  typedef typename LabelObjectType::AttributeType AttributeType;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ShapeRelabelLabelMapFilter, InPlaceLabelMapFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
/*  itkConceptMacro(InputEqualityComparableCheck,
    (Concept::EqualityComparable<InputImagePixelType>));
  itkConceptMacro(IntConvertibleToInputCheck,
    (Concept::Convertible<int, InputImagePixelType>));
  itkConceptMacro(InputOStreamWritableCheck,
    (Concept::OStreamWritable<InputImagePixelType>));*/
/** End concept checking */
#endif

  /**
   * Set/Get the order of labeling of the objects. By default, the objects with
   * the highest attribute values are labeled first. Setting ReverseOrdering to true
   * causes the object with the smallest attributes to be labeled first.
   */
  itkSetMacro(ReverseOrdering, bool);
  itkGetConstReferenceMacro(ReverseOrdering, bool);
  itkBooleanMacro(ReverseOrdering);

  /**
   * Set/Get the attribute to use.
   * Default is "Size".
   */
  itkGetConstMacro(Attribute, AttributeType);
  itkSetMacro(Attribute, AttributeType);
  void SetAttribute(const std::string & s)
  {
    this->SetAttribute( LabelObjectType::GetAttributeFromName(s) );
  }

protected:
  ShapeRelabelLabelMapFilter();
  ~ShapeRelabelLabelMapFilter() {}

  void GenerateData();

  typedef typename Superclass::LabelObjectContainerType LabelObjectContainerType;

  template< class TAttributeAccessor >
  void TemplatedGenerateData(const TAttributeAccessor &)
  {
    // Allocate the output
    this->AllocateOutputs();

    ImageType *output = this->GetOutput();

    const LabelObjectContainerType & labelObjectContainer = output->GetLabelObjectContainer();
    typedef typename LabelObjectType::Pointer LabelObjectPointer;
    typedef std::vector< LabelObjectPointer > VectorType;

    ProgressReporter progress( this, 0, 2 * labelObjectContainer.size() );

    // Get the label objects in a vector, so they can be sorted
    VectorType labelObjects;
    labelObjects.reserve( labelObjectContainer.size() );
    for ( typename LabelObjectContainerType::const_iterator it = labelObjectContainer.begin();
          it != labelObjectContainer.end();
          it++ )
      {
      labelObjects.push_back(it->second);
      progress.CompletedPixel();
      }

    // Instantiate the comparator and sort the vector
    if ( m_ReverseOrdering )
      {
      std::sort( labelObjects.begin(), labelObjects.end(),
                 Functor::LabelObjectReverseComparator< LabelObjectType, TAttributeAccessor >() );
      }
    else
      {
      std::sort( labelObjects.begin(), labelObjects.end(),
                 Functor::LabelObjectComparator< LabelObjectType, TAttributeAccessor >() );
      }
    //   progress.CompletedPixel();

    // and put back the objects in the map
    output->ClearLabels();
    unsigned int label = 0;
    typename VectorType::const_iterator it2 = labelObjects.begin();
    while ( it2 != labelObjects.end() )
      {
      // Avoid the background label if it is used
      if ( label == output->GetBackgroundValue() )
        {
        label++;
        }
      ( *it2 )->SetLabel(label);
      output->AddLabelObject(*it2);

      // Go to the next label
      label++;
      progress.CompletedPixel();

      it2++;
      }
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

  bool m_ReverseOrdering;

  AttributeType m_Attribute;
private:
  ShapeRelabelLabelMapFilter(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented
};                                          // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShapeRelabelLabelMapFilter.txx"
#endif

#endif
