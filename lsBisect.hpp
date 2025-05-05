#pragma once

#include <vcLogger.hpp>

template <class FUNC> class lsBisect {
  FUNC function;
  double a;
  double b;
  double result;
  unsigned numIter = 0;
  double eps = 1e-6;
  unsigned maxIterations = 100;

public:
  lsBisect() {}

  lsBisect(FUNC &func, double lower = 0., double upper = 1.)
      : function(func), a(lower), b(upper) {}

  void setFunc(FUNC &func) { function = func; }

  void setLimits(double lower, double upper) {
    a = lower;
    b = upper;
  }

  /// Set the x-accuracy to which the algorithm should go to
  /// before finishing. Defaults to 1e-6
  void setCriterion(double epsilon) { eps = epsilon; }

  /// Maximum number of iterations before aborting to not create
  /// infinite loop. Defaults to 100
  void setMaxIterations(unsigned maximumIterations) {
    maxIterations = maximumIterations;
  }

  double getRoot() const { return result; }

  void apply() {
    double funcA = function(a);

    if (funcA * function(b) > 0) {
      viennacore::Logger::getInstance().addError(
          "There is no sign change between lower and upper limit!");
    }

    result = a;
    for (; numIter < maxIterations && (b - a) > eps; ++numIter) {
      // take middle
      result = (a + b) / 2.;

      if (function(result) * funcA < 0) {
        b = result;
      } else {
        a = result;
        funcA = function(a);
      }
    }
  }
};