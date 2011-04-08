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
#ifndef __itkBorderQuadEdgeMeshFilter_h
#define __itkBorderQuadEdgeMeshFilter_h

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"

namespace itk
{
/**
 * \class BorderQuadEdgeMeshFilter
 * \brief Transform one border of a QuadEdgeMesh into either a circle
 * (conformal) or a square (arclength-wise).
 *
 * This class is one important step when computing a planar parameterization
 * of one mesh.
 *
 * If the input mesh has several boundaries, one can choose
 * the one which would be transformed via the variable m_BorderPick.
 *
 * \li <tt>m_BorderPick == Self::LONGEST</tt> refers to the boundary
 * \f$ b \f$ which satisfies:
 * \f[ b = \arg \max_{b^k} \sum_{i=1}^{N^k} \left\| x_{i}^k - x_{i+1}^k \right\| \f]
 *
 * \li <tt>m_BorderPick == Self::LARGEST</tt> refers to the boundary
 * \f$ b \f$ which satisfies:
 * \f[ b = \arg \max_{b^k} N^k \f]
 *
 * \sa ParameterizationQuadEdgeMeshFilter
 */
template< class TInputMesh, class TOutputMesh=TInputMesh >
class ITK_EXPORT BorderQuadEdgeMeshFilter:
  public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  /** Basic types. */
  typedef BorderQuadEdgeMeshFilter    Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh,
                                            TOutputMesh >
                                      Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  typedef TInputMesh                                  InputMeshType;
  typedef typename InputMeshType::ConstPointer        InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType        InputCoordRepType;
  typedef typename InputMeshType::PointType           InputPointType;
  typedef typename InputMeshType::Traits              InputTraits;
  typedef typename InputMeshType::PointIdentifier     InputPointIdentifier;
  typedef typename InputMeshType::QEType              InputQEType;
  typedef typename InputQEType::IteratorGeom          InputIteratorGeom;
  typedef typename InputMeshType::VectorType          InputVectorType;
  typedef typename InputMeshType::EdgeListType        InputEdgeListType;
  typedef typename InputMeshType::EdgeListPointerType InputEdgeListPointerType;
  typedef typename InputEdgeListType::iterator        InputEdgeListIterator;
  typedef typename InputMeshType::EdgeCellType        InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType     InputPolygonCellType;
  typedef typename InputMeshType::PointIdList         InputPointIdList;
  typedef typename InputMeshType::PointsContainer     InputPointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator
  InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstIterator
  InputCellsContainerConstIterator;

  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::Traits          OutputTraits;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::QEType          OutputQEType;
  typedef typename OutputMeshType::VectorType      OutputVectorType;
  typedef typename OutputMeshType::EdgeListType    OutputEdgeListType;
  typedef typename OutputMeshType::EdgeCellType    OutputEdgeCellType;
  typedef typename OutputMeshType::PolygonCellType OutputPolygonCellType;
  typedef typename OutputMeshType::PointIdList     OutputPointIdList;
  typedef typename OutputMeshType::PointsContainer OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerConstIterator
  OutputPointsContainerConstIterator;
  typedef typename OutputMeshType::CellsContainerConstIterator
  OutputCellsContainerConstIterator;

  itkNewMacro(Self);
  itkTypeMacro(BorderQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter);
  itkStaticConstMacro(PointDimension, unsigned int,
                      InputTraits::PointDimension);

  typedef std::vector< InputPointType >                           InputVectorPointType;
  typedef std::map< InputPointIdentifier, OutputPointIdentifier > MapPointIdentifier;
  typedef typename MapPointIdentifier::iterator                   MapPointIdentifierIterator;

  typedef QuadEdgeMeshBoundaryEdgesMeshFunction< InputMeshType > BoundaryRepresentativeEdgesType;
  typedef typename BoundaryRepresentativeEdgesType::Pointer      BoundaryRepresentativeEdgesPointer;

public:

  enum BorderTransformType {
    SQUARE_BORDER_TRANSFORM = 0,
    DISK_BORDER_TRANSFORM
    };

  enum BorderPickType {
    LONGEST = 0,
    LARGEST
    };

  itkSetMacro(TransformType, BorderTransformType);
  itkGetConstMacro(TransformType, BorderTransformType);

  itkSetMacro( BorderPick, BorderPickType );
  itkGetConstMacro( BorderPick, BorderPickType );

  itkSetMacro(Radius, InputCoordRepType);
  itkGetConstMacro(Radius, InputCoordRepType);

  void ComputeTransform();

  MapPointIdentifier GetBoundaryPtMap();

  InputVectorPointType GetBorder();

protected:
  /** \brief Constructor */
  BorderQuadEdgeMeshFilter();

  /** \brief Destructor */
  ~BorderQuadEdgeMeshFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const;

  BorderTransformType m_TransformType;
  BorderPickType m_BorderPick;

  InputCoordRepType m_Radius;

  InputVectorPointType m_Border;

  MapPointIdentifier m_BoundaryPtMap;

  void GenerateData();

  void ComputeBoundary();

  InputQEType* ComputeLongestBorder();

  InputQEType* ComputeLargestBorder();

  void DiskTransform();

  InputPointType GetMeshBarycentre();

  InputCoordRepType RadiusMaxSquare();

  void ArcLengthSquareTransform();

private:
  /** Not implemented */
  BorderQuadEdgeMeshFilter(const Self &);

  /** Not implemented */
  void operator=(const Self &);
};
}
#include "itkBorderQuadEdgeMeshFilter.txx"

#endif