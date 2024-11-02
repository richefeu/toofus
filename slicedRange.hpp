// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef SLICEDRANGE_HPP
#define SLICEDRANGE_HPP

#include <cmath>

template <typename T> class slicedRange {
public:
  /**
   * @brief Default constructor.
   *
   * Initializes the range with a minimum value of 0, a maximum value of 0, and 0 slices.
   */
  slicedRange() {
    vmin = 0;
    n = 0;
    scale = 1.0;
  }

  /**
   * @brief Constructor.
   *
   * Initializes the range with a minimum value of Vmin, a maximum value of Vmax, and N slices.
   *
   * @param Vmin The minimum value of the range.
   * @param Vmax The maximum value of the range.
   * @param N The number of slices.
   */
  slicedRange(T Vmin, T Vmax, int N) : vmin(Vmin), n(N) { scale = (T)n / (Vmax - Vmin); }

  /**
   * @brief Sets the minimum value, maximum value, and number of slices for the range.
   *
   * @param Vmin The minimum value of the range.
   * @param Vmax The maximum value of the range.
   * @param N The number of slices.
   */
  void set(T Vmin, T Vmax, int N) {
    vmin = Vmin;
    n = N;
    scale = (T)n / (Vmax - Vmin);
  }

  /**
   * @brief Sets the minimum and maximum values and the number of slices for the range.
   *
   * This method is a convenience wrapper that calls the set method with the given parameters.
   *
   * @param Vmin The minimum value of the range.
   * @param Vmax The maximum value of the range.
   * @param N The number of slices.
   */
  void set_MinMaxNb(T Vmin, T Vmax, int N) { set(Vmin, Vmax, N); }

  /**
   * @brief Sets the minimum value, the width of the slice, and the number of slices for the range.
   *
   * This method is a convenience wrapper that calls the set method with the given parameters.
   *
   * @param Vmin The minimum value of the range.
   * @param W The width of each slice.
   * @param N The number of slices.
   */
  void set_MinWidthNb(T Vmin, T W, int N) { set(Vmin, Vmin + (double)N * W, N); }

  /**
   * @brief Returns the ID of a given value in the range.
   *
   * The ID is the index of the slice in which the value is located.
   *
   * @param value The value to be located in the range.
   *
   * @return The ID of the value in the range. Returns a negative value if not in range.
   */
  int getID(T value) {
    int i = (int)floor(scale * (value - vmin));
    if (i >= n) {
      i = -i;
    }
    return i;
  }

  /**
   * @brief Returns the step size of each slice in the range.
   *
   * This method calculates the width of each slice based on the scale,
   * which is the inverse of the bin width. The step size is effectively
   * the distance between consecutive slices.
   *
   * @return The step size of each slice.
   */
  T getStep() const { return (1.0f / scale); }

  /**
   * @brief Returns the total number of slices in the range.
   *
   * This method provides the number of slices that the range is divided into.
   *
   * @return The total number of slices in the range.
   */
  int getNumberOfSlices() const { return n; }

  T getLeftValue() const { return vmin; }
  /**
   * @brief Returns the leftmost value of the range.
   *
   * This method retrieves the minimum value of the range, which is the starting point
   * for the slices.
   *
   * @return The leftmost value of the range.
   */

  /**
   * @brief Returns the rightmost value of the range.
   *
   * This method retrieves the maximum value of the range, which is the ending point
   * for the slices.
   *
   * @return The rightmost value of the range.
   */
  T getRightValue() const { return vmin + (double)n / scale; }

protected:
  T vmin;  ///< left value
  int n;   ///< Number of slices
  T scale; ///< inverse of the bin width
};

#endif /* end of include guard: SLICEDRANGE_HPP */
