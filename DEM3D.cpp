#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include "BoschProcess.hpp"
#include "MakeMask.hpp"

using namespace viennals;

int main() {
  omp_set_num_threads(16);

  constexpr int D = 3;
  typedef double NumericType;
  double gridDelta = 0.125;

  double extent = 12;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = BoundaryConditionEnum::INFINITE_BOUNDARY;

  auto mask = SmartPointer<Domain<NumericType, D>>::New(bounds, boundaryCons,
                                                        gridDelta);

  auto levelSet = SmartPointer<Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  std::array<NumericType, 3> maskOrigin = {};
  NumericType maskRadius = extent / 2.0;

  MakeMask<NumericType, D> maskCreator(levelSet, mask);
  maskCreator.setMaskOrigin(maskOrigin);
  maskCreator.setMaskRadius(maskRadius);
  maskCreator.apply();

  std::cout << "Output initial" << std::endl;
  auto mesh = SmartPointer<Mesh<NumericType>>::New();

  //   ToMesh<NumericType, D>(levelSet, mesh).apply();
  //   VTKWriter(mesh, "Surface_i_p.vtp").apply();
  ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  VTKWriter(mesh, "Surface_i.vtp").apply();
  //   ToMesh<NumericType, D>(mask, mesh).apply();
  //   VTKWriter(mesh, "Surface_m_p.vtp").apply();
  ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
  VTKWriter(mesh, "Surface_m.vtp").apply();

  NumericType bottomFraction = 0.7;
  NumericType etchRate = -1.86;

  BoschProcess<NumericType, D> processKernel(levelSet, mask);
  processKernel.setNumCycles(19);
  processKernel.setIsotropicRate(etchRate * 1.15);
  processKernel.setCycleEtchDepth(etchRate);
  processKernel.setStartWidth(2 * maskRadius);
  processKernel.setBottomWidth(2 * maskRadius * bottomFraction);
  processKernel.setStartOfTapering(-10);
  processKernel.setLateralEtchRatio(0.5);

  auto start = std::chrono::high_resolution_clock::now();
  processKernel.apply();
  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Geometric advect took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  std::cout << "Final structure has " << levelSet->getNumberOfPoints()
            << " LS points" << std::endl;

  // levelSet->print();
  ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  VTKWriter(mesh, "surface.vtp").apply();
  // ToMesh<NumericType, D>(levelSet, mesh).apply();
  // VTKWriter(mesh, "points-1.vtp").apply();

  std::cout << "Making volume output...(this may take a while)" << std::endl;

  auto volumeMeshing =
      SmartPointer<WriteVisualizationMesh<NumericType, D>>::New();
  volumeMeshing->insertNextLevelSet(mask);
  volumeMeshing->insertNextLevelSet(levelSet);
  volumeMeshing->setFileName("bosch");
  volumeMeshing->apply();

  return 0;
}