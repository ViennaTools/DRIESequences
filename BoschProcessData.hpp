#pragma once

#include <array>
#include <limits>

template <class T> struct BoschProcessDataType {
  unsigned numCycles;
  T isoRate;
  T startWidth;
  T bottomWidth;
  T taperStart = std::numeric_limits<T>::max();
  T topOffset = 0;
  std::array<T, 3> maskOrigin = {};
  bool sidewallTapering = true;
  bool scallopDecrease = true;
  double depthPerCycle = 0;
  unsigned numTaperCycles = 0;
  double taperRatio = 0.;
  double trenchBottom = 0.;
  double gridDelta = 0.;
  unsigned sausageCycle = 0;
  double sausageEtchRate = 0.;
  bool isWallTapering = true;
  double lateralRatio = 0.0;
};