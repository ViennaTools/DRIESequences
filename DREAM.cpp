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

double reFromAt(double a_t) {
  static constexpr double p0 = 1.17506441;
  static constexpr double p1 = 0.61536308;
  static constexpr double p2 = -0.42438527;
  static constexpr double t0 = p1 / p0 - p2;
  static constexpr double tm = p1 / (p0 - 1) - p2;

  if (a_t <= t0) {
    return 0.;
  } else if (a_t >= tm) {
    return 1.;
  } else {
    return p0 - p1 / (p2 + a_t);
  }
}

using namespace viennals;

int main() {
  omp_set_num_threads(16);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 0.025; // 0.125;

  double extent = 3;
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
  NumericType maskRadius = 0.4;

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

  // std::vector<NumericType> bottomFractions{
  //     0.13948421397683908, 0.41944813667047265, 0.6400617331600338,
  //     0.8183890372615938, 0.909870665101959};

  // Take average from etch rate measurements in Fig3.11d from ChangThesis2018
  NumericType etchRate =
      -(46 + 42 + 44 * 2) / (4 * 119.); //-0.3697478991596639; // 0.39

  // Calculate bottom fractions from model, based on ash times; last one is just
  // maximum
  std::vector<NumericType> ashTimes{0., 1.0, 1.2, 1.5, 2.0, 2.5, 4.0};
  std::vector<NumericType> bottomFractions;
  for (auto it : ashTimes) {
    bottomFractions.push_back(reFromAt(it));
    std::cout << bottomFractions.back() << std::endl;
  }

  BoschProcess<NumericType, D> processKernel;
  processKernel.setMask(mask);
  processKernel.setNumCycles(100);
  processKernel.setIsotropicRate(etchRate *
                                 0.6); // * 1.2 / 2 because it is a radius
  processKernel.setCycleEtchDepth(etchRate);
  processKernel.setStartWidth(2 * maskRadius);
  processKernel.setLateralEtchRatio(0.5);

  for (auto it : bottomFractions) {
    auto substrate = SmartPointer<Domain<NumericType, D>>::New(levelSet);
    processKernel.setSubstrate(substrate);
    processKernel.setBottomWidth(2 * maskRadius * it);
    processKernel.setStartOfTapering(-24.5);

    std::cout << "r_e: " << it << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    processKernel.apply();
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Geometric advect took: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop -
                                                                       start)
                     .count()
              << " ms" << std::endl;
    std::cout << "Final structure has " << substrate->getNumberOfPoints()
              << " LS points" << std::endl;

    // levelSet->print();
    ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    std::ostringstream out;
    out.precision(2);
    out << std::fixed << it;
    VTKWriter(mesh, "surface" + out.str() + ".vtp").apply();
    // ToMesh<NumericType, D>(levelSet, mesh).apply();
    // VTKWriter(mesh, "points-1.vtp").apply();

    std::cout << "Making volume output..." << std::endl;

    auto volumeMeshing =
        SmartPointer<WriteVisualizationMesh<NumericType, D>>::New();
    volumeMeshing->insertNextLevelSet(mask);
    volumeMeshing->insertNextLevelSet(substrate);
    volumeMeshing->setFileName("bosch" + out.str());
    volumeMeshing->apply();
  }

  return 0;
}