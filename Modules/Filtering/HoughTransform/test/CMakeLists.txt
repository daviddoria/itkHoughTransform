itk_module_test()
set(ITK-HoughTransformTests
itkHoughTransform2DLinesImageFilterTest.cxx
)

CreateTestDriver(ITK-HoughTransform  "${ITK-HoughTransform-Test_LIBRARIES}" "${ITK-HoughTransformTests}")

add_test(NAME itkHoughTransform2DLinesImageFilterTest
      COMMAND ITK-HoughTransformTestDriver itkHoughTransform2DLinesImageFilterTest)
