#pragma once

#include <lsSmartPointer.hpp>

template <class T, int D> class MakeMask {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;
  LSPtrType substrate;
  LSPtrType mask;

  std::array<T, 3> maskOrigin = {};
  T maskRadius = 0;
  T maskHeight = 2.;

  unsigned numberOfHoles = 1;
  double scallopSpacing = 1.5;

public:
  MakeMask(LSPtrType passedSubstrate, LSPtrType passedMask)
      : substrate(passedSubstrate), mask(passedMask) {}

  void setMaskOrigin(std::array<T, 3> &origin) { maskOrigin = origin; }

  void setMaskRadius(T radius) { maskRadius = radius; }

  void apply() {
    auto &grid = substrate->getGrid();
    auto &boundaryCons = grid.getBoundaryConditions();
    auto gridDelta = grid.getGridDelta();
    T extent = grid.getGridExtent(0) * gridDelta;

    // create mask
    {
      T normal[3] = {0., (D == 2) ? 1. : 0., (D == 3) ? 1. : 0.};
      T origin[3] = {0., (D == 2) ? maskHeight : 0.,
                     (D == 3) ? maskHeight : 0.};

      lsMakeGeometry<T, D>(mask,
                           lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
          .apply();
      normal[D - 1] = -1.0;
      origin[D - 1] = 0.;
      auto maskBottom = lsSmartPointer<lsDomain<T, D>>::New(grid);
      lsMakeGeometry<T, D>(maskBottom,
                           lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
          .apply();
      lsBooleanOperation<T, D>(mask, maskBottom,
                               lsBooleanOperationEnum::INTERSECT)
          .apply();

      // auto mesh = lsSmartPointer<lsMesh>::New();
      // lsToMesh<T, D>(mask, mesh).apply();
      // lsVTKWriter(mesh, "Plane.vtk").apply();

      auto maskHole = lsSmartPointer<lsDomain<T, D>>::New(grid);

      if constexpr (D == 3) {
        maskOrigin[2] = origin[2] - gridDelta;
        // maskRadius = extent / 2.0;
        double axis[3] = {0.0, 0.0, 1.0};
        lsMakeGeometry<T, D>(maskHole,
                             lsSmartPointer<lsCylinder<T, D>>::New(
                                 maskOrigin.data(), axis,
                                 maskHeight + 2 * gridDelta, maskRadius))
            .apply();
      } else {
        // double minScalar = origin[D-1] - 2 * gridDelta;
        // maskRadius = extent / 2.0;
        double min[3] = {-maskRadius, -gridDelta};
        double max[3] = {maskRadius, maskHeight + 2 * gridDelta};
        lsMakeGeometry<T, D>(maskHole,
                             lsSmartPointer<lsBox<T, D>>::New(min, max))
            .apply();
      }

      lsBooleanOperation<T, D>(mask, maskHole,
                               lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();

      // lsToSurfaceMesh<T, D>(mask, mesh).apply();
      // lsVTKWriter(mesh, "Mask.vtk").apply();

      // make substrate
      lsBooleanOperation<T, D>(maskBottom, lsBooleanOperationEnum::INVERT)
          .apply();
      substrate->deepCopy(maskBottom);
      lsBooleanOperation<T, D>(substrate, mask, lsBooleanOperationEnum::UNION)
          .apply();
    }
  }
};