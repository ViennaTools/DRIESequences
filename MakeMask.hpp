#pragma once

#include <lsDomain.hpp>

template <class T, int D> class MakeMask {
  using LSPtrType = viennals::SmartPointer<viennals::Domain<T, D>>;
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

      viennals::MakeGeometry<T, D>(
          mask,
          viennals::SmartPointer<viennals::Plane<T, D>>::New(origin, normal))
          .apply();
      normal[D - 1] = -1.0;
      origin[D - 1] = 0.;
      auto maskBottom =
          viennals::SmartPointer<viennals::Domain<T, D>>::New(grid);
      viennals::MakeGeometry<T, D>(
          maskBottom,
          viennals::SmartPointer<viennals::Plane<T, D>>::New(origin, normal))
          .apply();
      viennals::BooleanOperation<T, D>(
          mask, maskBottom, viennals::BooleanOperationEnum::INTERSECT)
          .apply();

      // auto mesh = viennals::SmartPointer<lsMesh<NumericType>>::New();
      // lsToMesh<T, D>(mask, mesh).apply();
      // lsVTKWriter(mesh, "Plane.vtp").apply();

      auto maskHole = viennals::SmartPointer<viennals::Domain<T, D>>::New(grid);

      if constexpr (D == 3) {
        maskOrigin[2] = origin[2] - gridDelta;
        // maskRadius = extent / 2.0;
        double axis[3] = {0.0, 0.0, 1.0};
        viennals::MakeGeometry<T, D>(
            maskHole, viennals::SmartPointer<viennals::Cylinder<T, D>>::New(
                          maskOrigin.data(), axis, maskHeight + 2 * gridDelta,
                          maskRadius))
            .apply();
      } else {
        // double minScalar = origin[D-1] - 2 * gridDelta;
        // maskRadius = extent / 2.0;
        double min[3] = {-maskRadius, -gridDelta};
        double max[3] = {maskRadius, maskHeight + 2 * gridDelta};
        viennals::MakeGeometry<T, D>(
            maskHole,
            viennals::SmartPointer<viennals::Box<T, D>>::New(min, max))
            .apply();
      }

      viennals::BooleanOperation<T, D>(
          mask, maskHole, viennals::BooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();

      // lsToSurfaceMesh<T, D>(mask, mesh).apply();
      // lsVTKWriter(mesh, "Mask.vtp").apply();

      // make substrate
      viennals::BooleanOperation<T, D>(maskBottom,
                                       viennals::BooleanOperationEnum::INVERT)
          .apply();
      substrate->deepCopy(maskBottom);
      viennals::BooleanOperation<T, D>(substrate, mask,
                                       viennals::BooleanOperationEnum::UNION)
          .apply();
    }
  }
};