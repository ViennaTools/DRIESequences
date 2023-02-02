#pragma once

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsSmartPointer.hpp>

template <class T, int D> class PillarMask {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;
  LSPtrType substrate;
  LSPtrType mask;

  std::array<T, 3> maskOrigin = {};
  T maskRadius = 0;
  T maskHeight = 1.;
  T lineDistance = 0;

public:
  PillarMask(LSPtrType passedSubstrate, LSPtrType passedMask)
      : substrate(passedSubstrate), mask(passedMask) {}

  void setMaskOrigin(std::array<T, 3> &origin) { maskOrigin = origin; }

  void setMaskRadius(T radius) { maskRadius = radius; }

  void setLineDistance(T distance) { lineDistance = distance; }

  void apply() {
    auto &grid = substrate->getGrid();
    auto &boundaryCons = grid.getBoundaryConditions();
    auto gridDelta = grid.getGridDelta();
    hrleVectorType<T, 6> domainBounds;
    for (unsigned i = 0; i < 3; ++i) {
      domainBounds[2 * i] = grid.getMinLocalCoordinate(i);
      domainBounds[2 * i + 1] = grid.getMaxLocalCoordinate(i);
    }

    // std::cout << "domainBounds: " << domainBounds << std::endl;
    // std::cout << "linDistance:  " << lineDistance << std::endl;
    // std::cout << "maskRadius:   " << maskRadius << std::endl;
    std::array<T, 3> maskVec;
    for (unsigned i = 0; i < 2; ++i) {
      maskVec[i] = domainBounds[2 * i] + 0.5 * lineDistance + maskRadius;
    }
    maskVec[2] = maskOrigin[2] - gridDelta;

    // while y is within bounds
    unsigned yLines = 1;
    double axis[3] = {0.0, 0.0, 1.0};
    while (maskVec[1] < domainBounds[3]) {
      // std::cout << "maskVec: " << maskVec[0] << ", " << maskVec[1] <<
      // std::endl; make single cylinder at maskVec
      auto maskSpot = lsSmartPointer<lsDomain<T, D>>::New(grid);
      lsMakeGeometry<T, 3>(maskSpot, lsSmartPointer<lsSphere<T, D>>::New(
                                         maskVec.data(), maskRadius))
          .apply();
      // lsMakeGeometry<T, 3>(maskSpot,
      //                       lsSmartPointer<lsCylinder<T, D>>::New(
      //                           maskVec.data(), axis,
      //                           maskHeight + 2 * gridDelta, maskRadius))
      //     .apply();
      // add cylinder to whole mask
      lsBooleanOperation<T, 3>(mask, maskSpot, lsBooleanOperationEnum::UNION)
          .apply();

      // advance maskVec in x direction
      maskVec[0] += 2 * (lineDistance + 2 * maskRadius);
      // if outside of x-range, advance y direction
      if (maskVec[0] > domainBounds[1]) {
        maskVec[0] = domainBounds[2] + (((yLines++) % 2) ? 1.5 : 0.5) *
                                           (lineDistance + 2 * maskRadius);
        maskVec[1] += lineDistance + 2 * maskRadius;
      }
    }

    // now make substrate level set
    auto maskBottom = lsSmartPointer<lsDomain<T, 3>>::New(mask);
    lsMakeGeometry<T, 3>(
        maskBottom, lsSmartPointer<lsPlane<T, 3>>::New(maskOrigin.data(), axis))
        .apply();
    substrate->deepCopy(maskBottom);
    lsBooleanOperation<T, 3>(substrate, mask, lsBooleanOperationEnum::UNION)
        .apply();
    lsBooleanOperation<T, 3>(maskBottom).apply(); // invert
    lsBooleanOperation<T, 3>(mask, maskBottom,
                             lsBooleanOperationEnum::INTERSECT)
        .apply();
  }
};