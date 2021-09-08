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
template<typename PointType=Point2D, size_t NumVerticesInLeaf = 20>
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
    std::unique_ptr<KDTreeNode<PointType, NumVerticesInLeaf>> leaf_; ///< initial leaf node containing all the points.
};


template<typename PointType, size_t NumVerticesInLeaf>
template<typename IterablePointStorage>
inline KDTree<PointType, NumVerticesInLeaf>::KDTree(const IterablePointStorage& points)
{
    PointType min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    PointType max(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    PointType extent(0, 0);
    for (auto& point : points)
    {
        getX(min) = std::min(getX(min), getX(point));
        getY(min) = std::min(getY(min), getY(point));
        getX(max) = std::max(getX(max), getX(point));
        getY(max) = std::max(getY(max), getY(point));
    }
    extent.x() = (getX(max) - getX(min)) * 1.1;
    extent.y() = (getY(max) - getY(min)) * 1.1;
    leaf_ = std::make_unique<KDTreeNode<PointType, NumVerticesInLeaf>>(min,
                                                                       extent,
                                                                       NODESPLIT::X);
    for (const auto& point : points)
    {
        leaf_->emplace(point);
    }
}

template<typename PointType, size_t NumVerticesInLeaf>
template<typename T>
inline PointType KDTree<PointType, NumVerticesInLeaf>::getClosestPoint(const T& point) const
{
    assert(leaf_ && "getClosestPoint requires at least one point in the KDTree.");
    PointType closestPoint;
    leaf_->getClosestPoint(closestPoint, point);

    return closestPoint;
}

template<typename PointType, size_t NumVerticesInLeaf>
inline size_t KDTree<PointType, NumVerticesInLeaf>::size()
{
    if (!leaf_)
    { return 0; }
    return leaf_->size();
}

template<typename PointType, size_t NumVerticesInLeaf>
void KDTree<PointType, NumVerticesInLeaf>::emplace(PointType&& point)
{
    if (!leaf_)
    {
        leaf_ = std::make_unique<KDTreeNode<PointType, NumVerticesInLeaf>>();
    }
    while (!leaf_->isPointInCurrentDomain(point))
    {
        leaf_ = std::make_unique<KDTreeNode<PointType, NumVerticesInLeaf>>(std::move(leaf_),
                                                                           leaf_->createParentLeaf(point));
    }
    leaf_->emplace(std::forward<PointType>(point));
}

template<typename PointType, size_t NumVerticesInLeaf>
template<class... Args>
inline void KDTree<PointType, NumVerticesInLeaf>::emplace(Args&& ... args)
{
    emplace(PointType(args ...));
}
}

#endif // header guard