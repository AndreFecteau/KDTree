/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at https://mozilla.org/MPL/2.0/.
/// author: Andre Fecteau <andre.fecteau1@gmail.com>

#ifndef KDTREE_POINT2D_H
#define KDTREE_POINT2D_H

namespace KDTree {

/// \struct Point2D
/// \details A structure providing coordinates in 2 dimensions.
struct Point2D
{

    /// Default Constructor
    Point2D() = default;

    /// Constructor
    /// \param x value for the first dimension
    /// \param y value for the second dimension
    Point2D(double x, double y) : x_(x), y_(y) {}

    /// Getter for the first dimension
    /// \return const ref to the value.
    [[nodiscard]] const double& x() const { return x_; }

    /// Getter for the second dimension
    /// \return const ref to the value.
    [[nodiscard]] const double& y() const { return y_; }

    /// Getter for the first dimension
    /// \return ref to the value.
    [[nodiscard]] double& x() { return x_; }

    /// Getter for the second dimension
    /// \return ref to the value.
    [[nodiscard]] double& y() { return y_; }

    double x_; ///< first dimension value
    double y_; ///< second dimension value
};

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline double& getX(Point2D& point) { return point.x(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline const double& getX(const Point2D& point) { return point.x(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline double& getY(Point2D& point) { return point.y(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline const double& getY(const Point2D& point) { return point.y(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline double& getX(std::pair<Point2D, size_t>& point) { return point.first.x(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline const double& getX(const std::pair<Point2D, size_t>& point) { return point.first.x(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline double& getY(std::pair<Point2D, size_t>& point) { return point.first.y(); }

/// Getter function that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
inline const double& getY(const std::pair<Point2D, size_t>& point) { return point.first.y(); }
}

#endif // header guard