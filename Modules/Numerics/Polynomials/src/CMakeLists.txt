set(ITK-Polynomials_SRC
itkMultivariateLegendrePolynomial.cxx
)

add_library(ITK-Polynomials ${ITK-Polynomials_SRC})
target_link_libraries(ITK-Polynomials  ${ITK-Common_LIBRARIES} ${TEST_DEPENDS_LIBRARIES} ${ITK-TestKernel_LIBRARIES})
itk_module_target(ITK-Polynomials)
