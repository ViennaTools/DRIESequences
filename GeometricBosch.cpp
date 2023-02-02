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

int main() {
  omp_set_num_threads(16);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 0.25;

  double extent = 12;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
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
  NumericType maskRadius = extent / 4.0;

  MakeMask<NumericType, D> maskCreator(levelSet, mask);
  maskCreator.setMaskOrigin(maskOrigin);
  maskCreator.setMaskRadius(maskRadius);
  maskCreator.apply();

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

  NumericType bottomFraction = 0.2;
  NumericType etchRate = -2.0;

  BoschProcess<NumericType, D> processKernel(levelSet, mask);
  processKernel.setNumCycles(200);
  processKernel.setIsotropicRate(etchRate);
  processKernel.setStartWidth(2 * maskRadius);
  processKernel.setBottomWidth(2 * maskRadius * bottomFraction);
  processKernel.setStartOfTapering(-190);
  processKernel.setTopOffset(etchRate / 3);

  processKernel.apply();

  // levelSet->print();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "surface.vtp").apply();
  // lsToMesh<NumericType, D>(levelSet, mesh).apply();
  // lsVTKWriter(mesh, "points-1.vtp").apply();

  std::cout << "Making volume output..." << std::endl;

  auto volumeMeshing =
      lsSmartPointer<lsWriteVisualizationMesh<NumericType, D>>::New();
  volumeMeshing->insertNextLevelSet(mask);
  volumeMeshing->insertNextLevelSet(levelSet);
  volumeMeshing->setFileName("bosch");
  volumeMeshing->apply();

  return 0;
}