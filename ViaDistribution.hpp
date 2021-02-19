#pragma once

#include <lsGeometricAdvectDistributions.hpp>

#include "BoschProcessData.hpp"

template <class T, int D>
class ViaDistribution : public lsGeometricAdvectDistribution<T, D> {
  double getDepth(const std::array<hrleCoordType, 3> &initial) const {
    if (!isTapering || std::abs(data.taperStart) > std::abs(data.trenchBottom)) {
      return data.trenchBottom;
    }

    double radius = 0;
    for (unsigned i = 0; i < D - 1; ++i) {
      double dist = initial[i] - data.maskOrigin[i];
      radius += dist * dist;
    }
    radius = std::max(0.0, std::sqrt(radius) - data.bottomWidth);

    // adjust depth depending on radius from the middle of the via
    double depth =
        data.taperStart +
        std::max(1.0 - radius / (data.startWidth - data.bottomWidth), 0.0) *
            taperDepth;
    return depth;
  }

public:
  BoschProcessDataType<T> data;
  const double taperDepth;
  const bool isTapering;

  ViaDistribution(const BoschProcessDataType<T> &processData)
      : data(processData), taperDepth(data.trenchBottom - data.taperStart), isTapering(data.sidewallTapering) {}

  bool isInside(const std::array<hrleCoordType, 3> &initial,
                const std::array<hrleCoordType, 3> &candidate,
                double eps = 0.) const {
    for (unsigned i = 0; i < D - 1; ++i) {
      if (std::abs(candidate[i] - initial[i]) > (data.gridDelta + eps)) {
        return false;
      }
    }
    if (std::abs(candidate[D - 1] - initial[D - 1]) >
        std::abs(data.trenchBottom) + eps) {
      return false;
    }
    return true;
  }

  T getSignedDistance(const std::array<hrleCoordType, 3> &initial,
                      const std::array<hrleCoordType, 3> &candidate) const {
    T distance = std::numeric_limits<T>::lowest();
    for (unsigned i = 0; i < D - 1; ++i) {
      T vector = std::abs(candidate[i] - initial[i]);
      distance = std::max(vector - data.gridDelta, distance);
    }
    T vector = std::abs(candidate[D - 1] - initial[D - 1]);
    distance = std::max(vector - std::abs(getDepth(initial)), distance);

    return (data.trenchBottom < 0) ? -distance : distance;
  }

  std::array<hrleCoordType, 6> getBounds() const {
    std::array<hrleCoordType, 6> bounds = {};
    for (unsigned i = 0; i < D - 1; ++i) {
      bounds[2 * i] = -data.gridDelta * ((data.trenchBottom < 0) ? -1 : 1);
      bounds[2 * i + 1] = data.gridDelta * ((data.trenchBottom < 0) ? -1 : 1);
    }
    bounds[2 * (D - 1)] = -data.trenchBottom;
    bounds[2 * (D - 1) + 1] = data.trenchBottom;

    return bounds;
  }
};