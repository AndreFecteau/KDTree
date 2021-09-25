/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at https://mozilla.org/MPL/2.0/.
/// author: Andre Fecteau <andre.fecteau1@gmail.com>

#ifndef KDTREE_POINT2D_H
#define KDTREE_POINT2D_H

namespace KDTree {

/// \struct Point2D
/// \details A structure providing coordinates in 2 dimensions.
template<typename T>
struct Point2D
{

    /// Default Constructor
    Point2D() = default;

    /// Constructor
    /// \param x value for the first dimension
    /// \param y value for the second dimension
    Point2D(const T& x, const T& y) : x_(x), y_(y) {}

    /// Getter for the first dimension
    /// \return const ref to the value.
    [[nodiscard]] const T& x() const { return x_; }

    /// Getter for the second dimension
    /// \return const ref to the value.
    [[nodiscard]] const T& y() const { return y_; }

    /// Getter for the first dimension
    /// \return ref to the value.
    [[nodiscard]] T& x() { return x_; }

    /// Getter for the second dimension
    /// \return ref to the value.
    [[nodiscard]] T& y() { return y_; }

    T x_; ///< first dimension value
    T y_; ///< second dimension value
};

/// Getter struct that is over-loadable for different types without modifying the type.
/// \details Current implementation of the KDTree depends on this overload to get the first and second coordinates
///          of it's template parameter.
/// \tparam PointType The storage of CoordinateType
/// \tparam CoordinateType The type of the coordinates
/// \tparam index index of the coordinate
template<typename PointType, typename CoordinateType, size_t index>
struct Access;

template<typename PointType, typename CoordinateType>
struct Access<PointType, CoordinateType, 0>
{
    static void set(PointType& point, const CoordinateType& value) { point.x_ = value; }

    static const CoordinateType& get(const PointType& point) { return point.x_; }
};

template<typename PointType, typename CoordinateType>
struct Access<PointType, CoordinateType, 1>
{
    static void set(PointType& point, const CoordinateType& value) { point.y_ = value; }

    static const CoordinateType& get(const PointType& point) { return point.y_; }
};

template<typename CoordinateType>
struct Access<std::pair<Point2D<CoordinateType>, size_t>, CoordinateType, 0>
{
    static void
    set(std::pair<Point2D<CoordinateType>, size_t>& point, const CoordinateType& value) { point.first.x_ = value; }

    static const CoordinateType& get(const std::pair<Point2D<CoordinateType>, size_t>& point) { return point.first.x_; }
};

template<typename CoordinateType>
struct Access<std::pair<Point2D<CoordinateType>, size_t>, CoordinateType, 1>
{
    static void
    set(std::pair<Point2D<CoordinateType>, size_t>& point, const CoordinateType& value) { point.first.y_ = value; }

    static const CoordinateType& get(const std::pair<Point2D<CoordinateType>, size_t>& point) { return point.first.y_; }
};
}

#endif // header guard