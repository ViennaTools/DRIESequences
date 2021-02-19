#pragma once

#include <lsGeometricAdvectDistributions.hpp>

#include "BoschProcessData.hpp"

template <class T, int D>
class BoschDistribution : public lsGeometricAdvectDistribution<T, D> {
public:
  BoschProcessDataType<T> data;

  const T scallopTop = 0.; // change this to be adjustable
  const T numericEps = 1e-3;

  const T gradient;
  const T taperPerCycle;
  const T logDenom;
  const T deltaO2;
  const T zPrefactor;
  const T isoRate;

  double calcZ(double n) const {
    const double &x = data.taperRatio;
    const double frac = (1 - x) / (1 + x);
    return std::abs(data.depthPerCycle) / (1 + x) * (1 - std::pow(frac, n)) /
           (1 - frac);
  }

  double getRadius(double z) const {
    const T linearFactor = std::min(1 - gradient * (data.taperStart - z), 1.0);

    if (z > scallopTop)
      return 0.0;
    // if z is at the bottom of the trench, always use the maximum radius
    if (std::abs(z + data.gridDelta - data.trenchBottom) < deltaO2)
      return data.isoRate * linearFactor;

    // check if within isotropic cycle of sausage sequence
    if(data.sausageCycle > 0) {
      double zMod = z - scallopTop;
      zMod = std::fmod(std::abs(zMod), std::abs(data.sausageCycle * data.depthPerCycle));

      // if z is within gridDelta/2 of sausage cycle
      if(std::abs(zMod - std::abs((data.sausageCycle-1) * data.depthPerCycle)) < deltaO2 + numericEps) {
        // std::cout << "sausageIso: " << data.sausageEtchRate << std::endl;
        return data.sausageEtchRate;
      }
    }

    if (std::abs(z) < std::abs(data.taperStart) - deltaO2) {
      z -= scallopTop;
      z = std::fmod(std::abs(z), std::abs(data.depthPerCycle));

      // if z is within gridDelta/2 of centre of cycle, use 1 otherwise 0
      if (std::abs(z - std::abs(data.depthPerCycle) / 2) <
          deltaO2 + numericEps) {
        return data.isoRate;
      } else {
        return 0;
      }
    } else {
      z -= data.taperStart;
      z = std::abs(z);

      const double n_z = std::log(1 - (z * zPrefactor)) / logDenom;

      const double nearestZ = calcZ(std::round(n_z));

      double diff = std::abs(z - nearestZ);

      if (diff < deltaO2 + numericEps) {
        return data.isoRate * linearFactor;
      } else {
        return 0;
      }
    }
  }

  BoschDistribution(BoschProcessDataType<T> &processData)
      : data(processData),
        gradient((1.0 - data.bottomWidth / data.startWidth) /
                 std::abs(data.trenchBottom - data.taperStart)),
        taperPerCycle((1 - data.taperRatio) / (1 + data.taperRatio)),
        logDenom(std::log(taperPerCycle)), deltaO2(data.gridDelta / 2.),
        zPrefactor(std::abs(2 * data.taperRatio / data.depthPerCycle)),
        isoRate((data.sausageCycle > 0)?data.sausageEtchRate:data.isoRate) {}

  bool isInside(const std::array<hrleCoordType, 3> &initial,
                const std::array<hrleCoordType, 3> &candidate,
                double eps = 0.) const {
    hrleCoordType dot = 0.;
    for (unsigned i = 0; i < D; ++i) {
      double tmp = candidate[i] - initial[i];
      dot += tmp * tmp;
    }

    if (std::sqrt(dot) <= std::abs(isoRate) + eps)
      return true;
    else
      return false;
  }

  T getSignedDistance(const std::array<hrleCoordType, 3> &initial,
                      const std::array<hrleCoordType, 3> &candidate) const {
    T currentRadius = getRadius(initial[D - 1]);
    T currentRadius2 = currentRadius * currentRadius;

    std::array<hrleCoordType, 3> v = {};
    for (unsigned i = 0; i < D; ++i) {
      v[i] = std::abs(candidate[i] - initial[i]);
      // subtract half of the length in x,y to generate "lens" distribution
      if (i < D - 1)
        v[i] += data.lateralRatio * currentRadius * ((data.isoRate < 0) ? -1 : 1);
    }

    if (std::abs(currentRadius) <= data.gridDelta) {
      T distance =
          std::max(std::max(std::abs(v[0]), std::abs(v[1])), std::abs(v[2])) -
          std::abs(currentRadius);
      return (currentRadius > 0) ? distance : -distance;
    }

    T distance = std::numeric_limits<T>::max();
    for (unsigned i = 0; i < D; ++i) {
      T y = (v[(i + 1) % D]);
      T z = 0;
      if constexpr (D == 3)
        z = (v[(i + 2) % D]);
      T x = currentRadius2 - y * y - z * z;
      if (x < 0.)
        continue;
      T dirRadius = v[i] - std::sqrt(x);
      if (std::abs(dirRadius) < std::abs(distance))
        distance = dirRadius;
    }
    return (data.isoRate > 0) ? distance : -distance;
  }

  std::array<hrleCoordType, 6> getBounds() const {
    std::array<hrleCoordType, 6> bounds{};
    for (unsigned i = 0; i < D - 1; ++i) {
      bounds[2 * i] = -isoRate;
      bounds[2 * i + 1] = isoRate;
    }
    bounds[2 * (D - 1)] = -isoRate;
    bounds[2 * D - 1] = isoRate;
    return bounds;
  }
};