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
#include "PillarMask.hpp"

int main() {
  omp_set_num_threads(32);

  constexpr int D = 3;
  typedef double NumericType;
  double gridDelta = 0.05; // 0.125;

  NumericType maskRadius = 1.25 / 2.;
  NumericType lineDistance = 0.51;
  NumericType unitCellLength = (2 * maskRadius + lineDistance);
  double extent = 2 * unitCellLength;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto mask = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  auto levelSet = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  std::array<NumericType, 3> maskOrigin = {};

  if constexpr (D == 2) {
    MakeMask<NumericType, D> maskCreator(levelSet, mask);
    maskCreator.setMaskOrigin(maskOrigin);
    maskCreator.setMaskRadius(maskRadius);
    maskCreator.apply();
  }

  if constexpr (D == 3) {
    PillarMask<NumericType, D> maskCreator(levelSet, mask);
    maskCreator.setMaskOrigin(maskOrigin);
    maskCreator.setMaskRadius(maskRadius);
    maskCreator.setLineDistance(lineDistance);
    maskCreator.apply();
  }

  std::cout << "Output initial" << std::endl;
  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();

  //   lsToMesh<NumericType, D>(levelSet, mesh).apply();
  //   lsVTKWriter(mesh, "Surface_i_p.vtp").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "Surface_i.vtp").apply();
  //   lsToMesh<NumericType, D>(mask, mesh).apply();
  //   lsVTKWriter(mesh, "Surface_m_p.vtp").apply();
  lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
  lsVTKWriter(mesh, "Surface_m.vtp").apply();

  // Take average from etch rate measurements in Fig6.c from Chang2018
  NumericType etchRate = -0.25;

  BoschProcess<NumericType, D> processKernel;
  processKernel.setMask(mask);
  processKernel.setNumCycles(80);
  processKernel.setIsotropicRate(etchRate *
                                 0.6); // * 1.2 / 2 because it is a radius
  processKernel.setCycleEtchDepth(etchRate);
  processKernel.setStartWidth(2 * maskRadius);

  processKernel.setSubstrate(levelSet);
  processKernel.setBottomWidth(2 * maskRadius);
  // processKernel.setStartOfTapering(-24.5);
  processKernel.setSausageCycling(10);
  processKernel.setSausageCycleDepth(2 * etchRate);
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
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "surface.vtp").apply();
  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "points-1.vtp").apply();

  std::cout << "Making volume output..." << std::endl;

  auto volumeMeshing =
      lsSmartPointer<lsWriteVisualizationMesh<NumericType, D>>::New();
  volumeMeshing->insertNextLevelSet(mask);
  volumeMeshing->insertNextLevelSet(levelSet);
  volumeMeshing->setFileName("bosch");
  volumeMeshing->apply();

  return 0;
}