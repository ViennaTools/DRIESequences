#pragma once

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>

template <class T, int D> class PillarMask {
  using LSPtrType = viennals::SmartPointer<viennals::Domain<T, D>>;
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
    viennacore::VectorType<T, 6> domainBounds;
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
      auto maskSpot = LSPtrType::New(grid);
      viennals::MakeGeometry<T, 3>(
          maskSpot, viennals::SmartPointer<viennals::Sphere<T, D>>::New(
                        maskVec.data(), maskRadius))
          .apply();
      // viennals::MakeGeometry<T, 3>(maskSpot,
      //                       viennals::SmartPointer<lsCylinder<T, D>>::New(
      //                           maskVec.data(), axis,
      //                           maskHeight + 2 * gridDelta, maskRadius))
      //     .apply();
      // add cylinder to whole mask
      viennals::BooleanOperation<T, 3>(mask, maskSpot,
                                       viennals::BooleanOperationEnum::UNION)
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
    auto maskBottom = viennals::SmartPointer<viennals::Domain<T, 3>>::New(mask);
    viennals::MakeGeometry<T, 3>(
        maskBottom, viennals::SmartPointer<viennals::Plane<T, 3>>::New(
                        maskOrigin.data(), axis))
        .apply();
    substrate->deepCopy(maskBottom);
    viennals::BooleanOperation<T, 3>(substrate, mask,
                                     viennals::BooleanOperationEnum::UNION)
        .apply();
    viennals::BooleanOperation<T, 3>(maskBottom).apply(); // invert
    viennals::BooleanOperation<T, 3>(mask, maskBottom,
                                     viennals::BooleanOperationEnum::INTERSECT)
        .apply();
  }
};