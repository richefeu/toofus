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

#ifndef SVGPLOT_HPP
#define SVGPLOT_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

struct SVGPlot {
  double padding{50.0};               //< Padding between the plot area and the SVG canvas border
  double width{500.0};                //< Width of the SVG canvas
  double height{300.0};               //< Height of the SVG canvas
  std::string xlabel{"x label"};      //< Label for the x-axis
  std::string ylabel{"y label"};      //< Label for the y-axis
  double xRangeMin{0.0};              //< Minimum value of the x-axis range
  double xRangeMax{10.0};             //< Maximum value of the x-axis range
  double yRangeMin{-1.0};             //< Minimum value of the y-axis range
  double yRangeMax{1.0};              //< Maximum value of the y-axis range
  double xTickInterval{1.0};          //< Interval between the x-axis ticks
  double yTickInterval{1.0};          //< Interval between the y-axis ticks
  bool showGrid{true};                //< Flag to show or hide the grid lines
  bool showBox{true};                 //< Flag to show or hide the surrounding box
  bool showTickLabels{true};          //< Flag to show or hide the tick labels
  bool showTickMarks{true};           //< Flag to show or hide the tick marks
  bool dashedLine{false};             //< Flag to enable or disable dashed lines for the plot
  bool hasMarkers{true};              //< Flag to enable or disable markers for the plot
  std::string fillMarkers{"none"};    //< Fill color for the plot markers
  std::string strokeMarkers{"black"}; //< Stroke color for the plot markers
  double radiusMarkers{3};            //< Radius of the markers
  int stepMarkers{1};                 //< Step between the data points to show markers
  double axisLabelDistance{30};       //< Distance between the axis labels and the axis
  double tickLabelDistance{5};        //< Distance between the tick labels and the axis
  int xTicksPrecision{2};             //< Precision (number of digits) for the x-tick labels
  int yTicksPrecision{2};             //< Precision (number of digits) for the y-tick labels

  /**
   * @brief Adjusts the range and interval of the given xmin and xmax values to create a specified number of intervals.
   *
   * @param xmin The minimum value of the range.
   * @param xmax The maximum value of the range.
   * @param desiredNumIntervals The desired number of intervals in the adjusted range.
   *
   * This function adjusts the range and interval of the given xmin and xmax values to create a specified number of
   * intervals. It first calculates the interval based on the desired number of intervals and adjusts xmin and xmax to
   * the nearest interval boundaries. Then, it adjusts the interval to a nice round value and updates xmin and xmax with
   * the adjusted interval.
   *
   * If the desiredNumIntervals is less than or equal to 0, or if xmin is greater than or equal to xmax, the function
   * will print an error message and return without modifying the input values.
   */
  void adjustRangeAndInterval(double &xmin, double &xmax, int desiredNumIntervals) {
    if (desiredNumIntervals <= 0) {
      std::cerr << "Error: desiredNumIntervals must be greater than 0." << std::endl;
      return;
    }

    if (xmin >= xmax) {
      std::cerr << "Error: xmin must be less than xmax." << std::endl;
      return;
    }

    double rangeWidth = xmax - xmin;
    double interval = rangeWidth / desiredNumIntervals;

    // Calculate the new xmin and xmax values based on the desired number of intervals
    xmin = std::floor(xmin / interval) * interval;
    xmax = std::ceil(xmax / interval) * interval;

    // Update the range width and interval values
    rangeWidth = xmax - xmin;
    interval = rangeWidth / desiredNumIntervals;

    // Adjust the interval to a nice round value
    double power = std::floor(std::log10(interval));
    double rangeFactor = std::pow(10, power);
    interval = std::ceil(interval / rangeFactor) * rangeFactor;

    // Update xmin and xmax with the adjusted interval
    xmin = std::floor(xmin / interval) * interval;
    xmax = xmin + desiredNumIntervals * interval;
  }

  /**
   * @brief Adjusts the ranges and intervals for the x and y axes based on the given data vectors.
   *
   * @param xvec The data vector for the x axis.
   * @param yvec The data vector for the y axis.
   * @param nx The desired number of intervals for the x axis (default: 10).
   * @param ny The desired number of intervals for the y axis (default: 10).
   *
   * This function adjusts the ranges and intervals for the x and y axes based on the given data vectors. It first finds
   * the minimum and maximum values of the x and y data vectors. Then, it adjusts the ranges and intervals for both axes
   * using the adjustRangeAndInterval() method with the desired number of intervals. Finally, it updates the xRangeMin,
   * xRangeMax, yRangeMin, yRangeMax, xTickInterval, and yTickInterval variables with the adjusted values.
   *
   * If either of the input data vectors is empty, the function will print an error message and return without modifying
   * the range and interval values.
   */
  void adjustRanges(const std::vector<double> &xvec, const std::vector<double> &yvec, int nx = 10, int ny = 10) {
    if (xvec.empty() || yvec.empty()) {
      std::cerr << "Error: Input data vectors cannot be empty." << std::endl;
      return;
    }

    // Find the minimum and maximum values of x and y
    double minX = *std::min_element(xvec.begin(), xvec.end());
    double maxX = *std::max_element(xvec.begin(), xvec.end());
    double minY = *std::min_element(yvec.begin(), yvec.end());
    double maxY = *std::max_element(yvec.begin(), yvec.end());

    adjustRangeAndInterval(minX, maxX, nx);
    adjustRangeAndInterval(minY, maxY, ny);

    xRangeMin = minX;
    xRangeMax = maxX;
    yRangeMin = minY;
    yRangeMax = maxY;
    xTickInterval = (maxX - minX) / (double)nx;
    yTickInterval = (maxY - minY) / (double)ny;
  }

  void generate(const std::vector<double> &xvec, const std::vector<double> &yvec, const std::string &filename) {
    std::ofstream file(filename);

    if (!file) {
      std::cout << "Error opening file: " << filename << std::endl;
      return;
    }

    // Set up the SVG file
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width << "\" height=\"" << height
         << "\">" << std::endl;

    // Calculate the plot dimensions
    double xScale = (width - 2 * padding) / (xRangeMax - xRangeMin);
    double yScale = (height - 2 * padding) / (yRangeMax - yRangeMin);

    // Draw the grid
    if (showGrid) {
      for (double x = xRangeMin + xTickInterval; x < xRangeMax; x += xTickInterval) {
        double xPos = padding + (x - xRangeMin) * xScale;
        file << "<line x1=\"" << xPos << "\" y1=\"" << padding << "\" x2=\"" << xPos << "\" y2=\"" << height - padding
             << "\" stroke=\"lightgray\" stroke-dasharray=\"2,2\"/>" << std::endl;
      }

      for (double y = yRangeMin + yTickInterval; y < yRangeMax; y += yTickInterval) {
        double yPos = height - padding - (y - yRangeMin) * yScale;
        file << "<line x1=\"" << padding << "\" y1=\"" << yPos << "\" x2=\"" << width - padding << "\" y2=\"" << yPos
             << "\" stroke=\"lightgray\" stroke-dasharray=\"2,2\"/>" << std::endl;
      }
    }

    if (showBox) {
      file << "<rect x=\"" << padding << "\" y=\"" << padding << "\" width=\"" << width - 2 * padding << "\" height=\""
           << height - 2 * padding << "\" stroke=\"black\" fill=\"none\" />" << std::endl;
    }

    // Draw the x-axis
    file << "<line x1=\"" << padding << "\" y1=\"" << height - padding << "\" x2=\"" << width - padding << "\" y2=\""
         << height - padding << "\" stroke=\"black\"/>" << std::endl;

    // Draw the y-axis
    file << "<line x1=\"" << padding << "\" y1=\"" << padding << "\" x2=\"" << padding << "\" y2=\"" << height - padding
         << "\" stroke=\"black\"/>" << std::endl;

    // Draw the x-axis label
    file << "<text x=\"" << width / 2.0 << "\" y=\"" << height - padding + axisLabelDistance
         << "\" font-size=\"12px\" text-anchor=\"middle\">" << xlabel << "</text>" << std::endl;

    // Draw the y-axis label
    file << "<text x=\"" << padding / 2.0 << "\" y=\"" << height / 2.0
         << "\" font-size=\"12px\" text-anchor=\"middle\" transform=\"rotate(-90, " << padding - axisLabelDistance
         << ", " << height / 2.0 << ")\">" << ylabel << "</text>" << std::endl;

    // Draw the x-axis tick labels
    if (showTickLabels) {
      for (double x = xRangeMin; x <= xRangeMax; x += xTickInterval) {
        double xPos = padding + (x - xRangeMin) * xScale;
        double yPos = height - padding + tickLabelDistance + 10;
        std::stringstream xLabel;
        xLabel << std::fixed << std::setprecision(xTicksPrecision) << x;
        file << "<text x=\"" << xPos << "\" y=\"" << yPos << "\" font-size=\"10px\" text-anchor=\"middle\">"
             << xLabel.str() << "</text>" << std::endl;
      }
      if (showTickMarks) {
        for (double x = xRangeMin; x <= xRangeMax; x += xTickInterval) {
          double xPos = padding + (x - xRangeMin) * xScale;
          double yPos = height - padding;
          file << "<line x1=\"" << xPos << "\" y1=\"" << yPos << "\" x2=\"" << xPos << "\" y2=\"" << yPos + 5
               << "\" stroke=\"black\"/>" << std::endl;
        }
      }
    }

    // Draw the y-axis tick labels
    if (showTickLabels) {
      for (double y = yRangeMin; y <= yRangeMax; y += yTickInterval) {
        double xPos = padding - tickLabelDistance;
        double yPos = height - padding - (y - yRangeMin) * yScale;
        std::stringstream yLabel;
        yLabel << std::fixed << std::setprecision(yTicksPrecision) << y;
        file << "<text x=\"" << xPos << "\" y=\"" << yPos << "\" font-size=\"10px\" text-anchor=\"end\">"
             << yLabel.str() << "</text>" << std::endl;
      }
      if (showTickMarks) {
        for (double y = yRangeMin; y <= yRangeMax; y += yTickInterval) {
          double xPos = padding;
          double yPos = height - padding - (y - yRangeMin) * yScale;
          file << "<line x1=\"" << xPos - 5 << "\" y1=\"" << yPos << "\" x2=\"" << xPos << "\" y2=\"" << yPos
               << "\" stroke=\"black\"/>" << std::endl;
        }
      }
    }

    // Draw the plot line
    file << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2";
    if (dashedLine) {
      file << "\" stroke-dasharray=\"5,5";
    }
    file << "\" points=\"";
    for (size_t i = 0; i < xvec.size(); ++i) {
      double x = padding + (xvec[i] - xRangeMin) * xScale;
      double y = height - padding - (yvec[i] - yRangeMin) * yScale;
      file << x << "," << y << " ";
    }
    file << "\"/>" << std::endl;

    // Draw the plot markers
    if (hasMarkers) {
      for (size_t i = 0; i < xvec.size(); i += stepMarkers) {
        double x = padding + (xvec[i] - xRangeMin) * xScale;
        double y = height - padding - (yvec[i] - yRangeMin) * yScale;
        file << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << radiusMarkers << "\" stroke=\"" << strokeMarkers
             << "\" fill=\"" << fillMarkers << "\"/>" << std::endl;
      }
    }

    // Close the SVG file
    file << "</svg>" << std::endl;

    file.close();

    std::cout << "SVG plot generated successfully: " << filename << std::endl;
  }
};

#endif /* end of include guard: SVGPLOT_HPP */

#if 0

int main(int argc, char const *argv[]) {

  std::vector<double> xvec;
  std::vector<double> yvec;

  for (double x = -2 * M_PI; x < M_PI; x += M_PI / 100.0) {
    xvec.push_back(x);
    yvec.push_back(x * sin(x) + cos(x));
  }

  SVGPlot plt;
  plt.xlabel = "time [s]";
  plt.ylabel = "perf [-]";
  plt.showGrid = true;
  plt.showTickLabels = true;
  plt.dashedLine = true;
  plt.strokeMarkers = "blue";
  plt.fillMarkers = "none";
  plt.stepMarkers = 20;
  plt.radiusMarkers = 5.0;
  plt.adjustRanges(xvec, yvec, 5, 4);
  plt.generate(xvec, yvec, "plot.svg");

  plt.dashedLine = false;
  plt.hasMarkers = false;
  plt.xRangeMin = -6.4;
  plt.xRangeMax = 6.4;
  plt.xTickInterval = 0.8;
  plt.yRangeMin = -5.0;
  plt.yRangeMax = 2.0;
  plt.yTickInterval = 1.0;
  plt.yTicksPrecision = 0;
  plt.xTicksPrecision = 1;
  plt.padding = 45.0;
  plt.showGrid = false;
  plt.generate(xvec, yvec, "plot1.svg");
  return 0;
}

#endif