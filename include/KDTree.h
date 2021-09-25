/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at https://mozilla.org/MPL/2.0/.
/// author: Andre Fecteau <andre.fecteau1@gmail.com>

#ifndef KDTREE_KDTREE_H
#define KDTREE_KDTREE_H

#include <optional>
#include "Point2D.h"
#include "KDTreeNode.h"

namespace KDTree {

/// \class KDTree
/// \brief Simple tree structure with alternating half splitting leafs ... for now.
/// \details Simple tree structure
///          - Tree to incrementally add points to the structure.
///          - Get the closest point to a given input.
///          - Does not check for duplicates, expect unique points.
/// \tparam PointType the point type that the container will use.
/// \tparam NumVerticesInLeaf The number of points per leaf.
template<typename PointType = Point2D<double>, typename CoordinateType = double, size_t NumVerticesInLeaf = 20>
class KDTree
{
public:

    /// \brief Default constructor
    KDTree() = default;

    /// \brief Custom constructor
    /// \tparam IterablePointStorage An iterative storage containing points
    /// \param points storage containing points
    template<typename IterablePointStorage>
    explicit KDTree(const IterablePointStorage& points);

    /// \brief get the closest points to the input.
    /// \tparam T allows to get the closest point given any point representation having .x() and .y() functions.
    /// \param point The point you want to get the closest point from.
    template<typename T>
    [[nodiscard]] PointType getClosestPoint(const T& point) const;

    /// \brief get the number of points in the container.
    [[nodiscard]] size_t size();

    /// \brief Check if the container is empty.
    [[nodiscard]] bool empty() { return !leaf_; }

    /// \brief Insert a point in the storage.
    /// \param point The point to insert.
    void emplace(PointType&& point);

    /// \brief Insert a point in the storage.
    /// \param args Constructor arguments for a point.
    template<class... Args>
    void emplace(Args&& ... args);

private:
    std::unique_ptr<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>> leaf_; ///< initial leaf node containing all the points.
};


template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<typename IterablePointStorage>
inline KDTree<PointType, CoordinateType,
        NumVerticesInLeaf>::KDTree(const IterablePointStorage& points)
{
    PointType min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    PointType max(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    PointType extent(0, 0);
    for (auto& point : points)
    {
        Access<PointType, CoordinateType, 0>::get(min) = std::min(Access<PointType, CoordinateType, 0>::get(min),
                                                                  Access<PointType, CoordinateType, 0>::get(point));
        Access<PointType, CoordinateType, 1>::get(min) = std::min(Access<PointType, CoordinateType, 1>::get(min),
                                                                  Access<PointType, CoordinateType, 1>::get(point));
        Access<PointType, CoordinateType, 0>::get(max) = std::max(Access<PointType, CoordinateType, 0>::get(max),
                                                                  Access<PointType, CoordinateType, 0>::get(point));
        Access<PointType, CoordinateType, 1>::get(max) = std::max(Access<PointType, CoordinateType, 1>::get(max),
                                                                  Access<PointType, CoordinateType, 1>::get(point));
    }
    extent.x() =
            (Access<PointType, CoordinateType, 0>::get(max) - Access<PointType, CoordinateType, 0>::get(min)) * 1.1;
    extent.y() =
            (Access<PointType, CoordinateType, 1>::get(max) - Access<PointType, CoordinateType, 1>::get(min)) * 1.1;
    leaf_ = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(min,
                                                                                       extent,
                                                                                       NODESPLIT::X);
    for (const auto& point : points)
    {
        leaf_->emplace(point);
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<typename T>
inline PointType KDTree<PointType, CoordinateType, NumVerticesInLeaf>::getClosestPoint(const T& point) const
{
    assert(leaf_ && "getClosestPoint requires at least one point in the KDTree.");
    PointType closestPoint;
    leaf_->getClosestPoint(closestPoint, point);

    return closestPoint;
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline size_t KDTree<PointType, CoordinateType, NumVerticesInLeaf>::size()
{
    if (!leaf_)
    { return 0; }
    return leaf_->size();
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
void KDTree<PointType, CoordinateType, NumVerticesInLeaf>::emplace(PointType&& point)
{
    if (!leaf_)
    {
        leaf_ = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>();
    }
    while (!leaf_->isPointInCurrentDomain(point))
    {
        leaf_ = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(std::move(leaf_),
                                                                                           leaf_->createParentLeaf(
                                                                                                   point));
    }
    leaf_->emplace(std::forward<PointType>(point));
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<class... Args>
inline void KDTree<PointType, CoordinateType, NumVerticesInLeaf>::emplace(Args&& ... args)
{
    emplace(PointType(args ...));
}
}

#endif // header guard