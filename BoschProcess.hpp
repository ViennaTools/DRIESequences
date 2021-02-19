#include <array>

#include <lsDomain.hpp>
#include <lsSmartPointer.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include "BoschDistribution.hpp"
#include "BoschProcessData.hpp"
#include "ViaDistribution.hpp"
#include "lsBisect.hpp"

template <class T, int D> class BoschProcess {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;

  LSPtrType substrate;
  LSPtrType mask;

  BoschProcessDataType<T> processData;

  double taperRatioFromRe(double r_e) {
    auto lambda = [r_e, N_t = processData.numTaperCycles](double x) {
      return 1 - (1 - std::pow((1 - x) / (1 + x), N_t)) * (1 + x) - r_e;
    };

    lsBisect bisect(lambda, 1e-6, 1.0);
    bisect.setCriterion(1e-9);
    bisect.apply();
    return bisect.getRoot();
  }

  double zFromTaperRatio() {
    const double &x = processData.taperRatio;
    const double frac = (1 - x) / (1 + x);
    return processData.depthPerCycle / (1 + x) *
           (1 - std::pow(frac, processData.numTaperCycles-1)) / (1 - frac);
  }

public:
  BoschProcess() {}

  BoschProcess(LSPtrType passedSubstrate, LSPtrType passedMask)
      : substrate(passedSubstrate), mask(passedMask) {
    processData.gridDelta = substrate->getGrid().getGridDelta();
  }

  void setSubstrate(LSPtrType levelSet) {
    substrate = levelSet;
    processData.gridDelta = substrate->getGrid().getGridDelta();
  }

  void setMask(LSPtrType levelSet) { mask = levelSet; }

  void setNumCycles(unsigned numberOfCycles) {
    processData.numCycles = numberOfCycles;
  }

  void setIsotropicRate(T isotropicRate) {
    processData.isoRate = isotropicRate;
  }

  void setStartWidth(T widthOfTrenchTop) {
    processData.startWidth = widthOfTrenchTop / 2.;
  }

  void setBottomWidth(T widthOfTrenchBottom) {
    processData.bottomWidth = widthOfTrenchBottom / 2.;
  }

  void setStartOfTapering(T startOfTapering) {
    processData.taperStart = startOfTapering;
  }

  void setTopOffset(T offsetAtTopScallop) {
    processData.topOffset = offsetAtTopScallop;
  }

  void setMaskOrigin(std::array<T, 3> centreOfMask) {
    processData.maskOrigin = centreOfMask;
  }

  void setTapering(bool isTapering) {
    processData.sidewallTapering = isTapering;
  }

  void setScallopDecrease(bool isScallopDecrease) {
    processData.scallopDecrease = isScallopDecrease;
  }

  void setCycleEtchDepth(double cycleEtchDepth) {
    processData.depthPerCycle = cycleEtchDepth;
  }

  void setSausageCycling(unsigned nthIsSausageCycle) {
    processData.sausageCycle = nthIsSausageCycle;
  }

  void setSausageCycleDepth(double sausageCycleEtchRate) {
    processData.sausageEtchRate = sausageCycleEtchRate;
  }

  void setSidewallTapering(bool isSidewallTapering) {
    processData.isWallTapering = isSidewallTapering;
  }

  void setLateralEtchRatio(double ratioLateral) {
    processData.lateralRatio = std::max(std::min(1.0 - ratioLateral, 1.0), 0.0);
  }

  void apply() {
    // calculate all required values
    const double r_e = processData.bottomWidth / processData.startWidth;
    processData.trenchBottom =
        processData.depthPerCycle * processData.numCycles +
        processData.topOffset - processData.gridDelta / 2.;

    if (std::abs(r_e - 1.0) < 1e-3) {
      processData.taperStart = std::numeric_limits<T>::lowest();
    }

    // if there is tapering
    if (std::abs(processData.taperStart) < std::abs(processData.trenchBottom)) {
      unsigned numStraightCycles =
          std::ceil(processData.taperStart / processData.depthPerCycle);
      processData.numTaperCycles = processData.numCycles - numStraightCycles;

      processData.taperStart = processData.depthPerCycle * numStraightCycles +
                               processData.topOffset +
                               processData.depthPerCycle / 2.;

      processData.taperRatio = taperRatioFromRe(r_e);

      processData.trenchBottom = processData.taperStart + zFromTaperRatio();
    }

    std::cout << "d_c: " << processData.depthPerCycle << std::endl;
    std::cout << "N_t: " << processData.numTaperCycles << std::endl;
    std::cout << "L_t: " << processData.taperStart << std::endl;
    std::cout << "r_e: " << r_e << std::endl;
    std::cout << "x:   " << processData.taperRatio << std::endl;
    std::cout << "L_b: " << processData.trenchBottom << std::endl;

    auto dist = lsSmartPointer<ViaDistribution<T, D>>::New(processData);

    lsGeometricAdvect<T, D>(substrate, dist, mask).apply();

#ifndef NDEBUG
    {
      auto mesh = lsSmartPointer<lsMesh>::New();
      // lsToMesh<T, D>(substrate, mesh).apply();
      // lsVTKWriter(mesh, "points-0.vtk").apply();
      lsToSurfaceMesh<T, D>(substrate, mesh).apply();
      lsVTKWriter(mesh, "DEBUG_BoschProcess_0.vtk").apply();
      auto writer = lsWriteVisualizationMesh<T, D>();
      writer.insertNextLevelSet(mask);
      writer.insertNextLevelSet(substrate);
      writer.setFileName("DEBUG_BoschProcess_v1");
      writer.apply();
      // std::cout << "Making scallops" << std::endl;
    }
#endif

    // Now make scallops on the sidewalls
    auto boschDist = lsSmartPointer<BoschDistribution<T, D>>::New(processData);

    // perform geometric advection
    lsGeometricAdvect<T, D> fastAdvectKernel(substrate, boschDist, mask);

    fastAdvectKernel.apply();

#ifndef NDEBUG
    {
      // substrate->print();
      auto mesh = lsSmartPointer<lsMesh>::New();
      lsToSurfaceMesh<T, D>(substrate, mesh).apply();
      lsVTKWriter(mesh, "DEBUG_BoschProcess_1.vtk").apply();
      // lsToMesh<T, D>(substrate, mesh).apply();
      // lsVTKWriter(mesh, "points-1.vtk").apply();
      auto writer = lsWriteVisualizationMesh<T, D>();
      writer.insertNextLevelSet(mask);
      writer.insertNextLevelSet(substrate);
      writer.setFileName("DEBUG_BoschProcess_v2");
      writer.apply();
    }
#endif
  }
};