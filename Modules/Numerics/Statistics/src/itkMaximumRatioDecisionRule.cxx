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
#include "itkMaximumRatioDecisionRule.h"

namespace itk
{
MaximumRatioDecisionRule::MaximumRatioDecisionRule()
{}

void MaximumRatioDecisionRule::SetAPriori(APrioriVectorType & values)
{
  m_NumberOfClasses = values.size();
  m_APrioriRatioMatrix.set_size( values.size(), values.size() );
  APrioriVectorSizeType i, j;
  double                APrioriRatio;
  for ( i = 0; i < m_NumberOfClasses; i++ )
    {
    for ( j = 0; j < m_NumberOfClasses; j++ )
      {
      if ( values[i] > 0 )
        {
        APrioriRatio = (double)values[j]
                       / (double)values[i];
        }
      else
        {
        APrioriRatio = NumericTraits< double >::max();
        }
      m_APrioriRatioMatrix.put(i, j, APrioriRatio);
      }
    }
}

unsigned int
MaximumRatioDecisionRule::Evaluate(const VectorType & discriminantScores) const
{
  unsigned int i, j;
  double       temp;

  for ( i = 0; i < m_NumberOfClasses; i++ )
    {
    j = 0;
    while ( j < m_NumberOfClasses )
      {
      if ( j != i )
        {
        if ( discriminantScores[j] != 0.0 )
          {
          temp = discriminantScores[i] / discriminantScores[j];
          }
        else
          {
          temp = NumericTraits< double >::max();
          }

        if ( temp < m_APrioriRatioMatrix.get(i, j) )
          {
          break;
          }
        }

      ++j;

      if ( j == m_NumberOfClasses )
        {
        return i;
        }
      }
    }

  return i;
}

unsigned int
MaximumRatioDecisionRule::Evaluate(const ArrayType & discriminantScores) const
{
  unsigned int i, j;
  double       temp;

  for ( i = 0; i < m_NumberOfClasses; i++ )
    {
    j = 0;
    while ( j < m_NumberOfClasses )
      {
      if ( j != i )
        {
        if ( discriminantScores[j] != 0.0 )
          {
          temp = discriminantScores[i] / discriminantScores[j];
          }
        else
          {
          temp = NumericTraits< double >::max();
          }

        if ( temp < m_APrioriRatioMatrix.get(i, j) )
          {
          break;
          }
        }

      ++j;

      if ( j == m_NumberOfClasses )
        {
        return i;
        }
      }
    }

  return i;
}
} // end of namespace
