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
#ifndef __itkRigid3DTransform_h
#define __itkRigid3DTransform_h

#include <iostream>
#include "itkMatrixOffsetTransformBase.h"
#include "itkVersor.h"

namespace itk
{
/** \class Rigid3DTransform
 * \brief Rigid3DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rotation and translation in 3D space.
 * The transform is specified as a rotation matrix around a arbitrary center
 * and is followed by a translation.
 *
 * The parameters for this transform can be set either using individual Set
 * methods or in serialized form using SetParameters() and SetFixedParameters().
 *
 * The serialization of the optimizable parameters is an array of 12 elements.
 * The first 9 parameters represents the rotation matrix in row-major order
 * (where the column index varies the fastest). The last 3 parameters defines
 * the translation in each dimension.
 *
 * The serialization of the fixed parameters is an array of 3 elements defining
 * the center of rotation in each dimension.
 *
 * \ingroup Transforms
 * \ingroup ITK-Transform
 */
template< class TScalarType = double >
// type for scalars (float or double)
class ITK_EXPORT Rigid3DTransform:
  public MatrixOffsetTransformBase< TScalarType, 3, 3 >
{
public:
  /** Standard class typedefs. */
  typedef Rigid3DTransform                               Self;
  typedef MatrixOffsetTransformBase< TScalarType, 3, 3 > Superclass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

#ifdef ITKV3_COMPATIBILITY
  /** Run-time type information (and related methods).   */
  itkNewMacro(Self);
#endif

  /** Run-time type information (and related methods). */
  itkTypeMacro(Rigid3DTransform, MatrixOffsetTransformBase);

  /** Dimension of the space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 12);

  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::ParametersValueType       ParametersValueType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::OutputVectorValueType     OutputVectorValueType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::MatrixValueType           MatrixValueType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::TranslationType           TranslationType;
  typedef typename Superclass::OffsetType                OffsetType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 12 parameters. The first 9 represents the rotation
   * matrix is row-major order and the last 3 represents the translation.
   *
   * \warning The rotation matrix must be orthogonal to within a specified tolerance,
   * else an exception is thrown.
   *
   * \sa Transform::SetParameters()
   * \sa Transform::SetFixedParameters() */
  virtual void SetParameters(const ParametersType & parameters);

  /** Directly set the rotation matrix of the transform.
   * \warning The input matrix must be orthogonal to within a specified tolerance,
   * else an exception is thrown.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix(const MatrixType & matrix);

  /**
   * Get rotation Matrix from an Rigid3DTransform
   *
   * This method returns the value of the rotation of the
   * Rigid3DTransform.
   *
   * \deprecated Use GetMatrix instead
   */
  const MatrixType & GetRotationMatrix()
  { return this->GetMatrix(); }

  /**
   * Set the rotation Matrix of a Rigid3D Transform
   *
   * This method sets the 3x3 matrix representing a rotation
   * in the transform.  The Matrix is expected to be orthogonal
   * with a certain tolerance.
   *
   * \deprecated Use SetMatrix instead
   *
   */
  virtual void SetRotationMatrix(const MatrixType & matrix)
  { this->SetMatrix(matrix); }

  /**
   * Compose the transformation with a translation
   *
   * This method modifies self to include a translation of the
   * origin.  The translation is precomposed with self if pre is
   * true, and postcomposed otherwise.
   */
  void Translate(const OffsetType & offset, bool pre = false);

#ifdef ITKV3_COMPATIBILITY
/** Get an inverse of this transform. */
  bool GetInverse(Self *inverse) const
  {
  return this->Superclass::GetInverse(inverse);
  }

/** Return an inverse of this transform. */
virtual InverseTransformBasePointer GetInverseTransform() const
  {
  Pointer inv = New();
  return this->GetInverse(inv) ? inv.GetPointer() : NULL;
  }
#endif

  /**
   * Back transform by an affine transformation
   *
   * This method finds the point or vector that maps to a given
   * point or vector under the affine transformation defined by
   * self.  If no such point exists, an exception is thrown.
   *
   * \deprecated Please use GetInverseTransform and then call the forward
   *   transform using the result.
   *
   */
  InputPointType      BackTransform(const OutputPointType
                                    & point) const;

  InputVectorType     BackTransform(const OutputVectorType
                                    & vector) const;

  InputVnlVectorType  BackTransform(const OutputVnlVectorType
                                    & vector) const;

  InputCovariantVectorType BackTransform(const OutputCovariantVectorType
                                         & vector) const;

  /**
   * Utility function to test if a matrix is orthogonal within a specified
   * tolerance
   */
  bool MatrixIsOrthogonal(const MatrixType & matrix, double tol = 1e-10);

protected:
  Rigid3DTransform(unsigned int spaceDim,
                   unsigned int paramDim);
  Rigid3DTransform(const MatrixType & matrix,
                   const OutputVectorType & offset);
  Rigid3DTransform();
  ~Rigid3DTransform();

  /**
   * Print contents of an Rigid3DTransform
   */
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  Rigid3DTransform(const Self &); //purposely not implemented
  void operator=(const Self &);   //purposely not implemented
};                                //class Rigid3DTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_Rigid3DTransform(_, EXPORT, TypeX, TypeY)                \
  namespace itk                                                               \
  {                                                                           \
  _( 1 ( class EXPORT Rigid3DTransform< ITK_TEMPLATE_1 TypeX > ) )            \
  namespace Templates                                                         \
  {                                                                           \
  typedef Rigid3DTransform< ITK_TEMPLATE_1 TypeX > Rigid3DTransform##TypeY; \
  }                                                                           \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkRigid3DTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkRigid3DTransform.txx"
#endif

#endif /* __itkRigid3DTransform_h */
