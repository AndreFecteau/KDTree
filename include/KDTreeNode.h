/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at https://mozilla.org/MPL/2.0/.
/// author: Andre Fecteau <andre.fecteau1@gmail.com>

#ifndef KDTREE_KDTTREENODE_H
#define KDTREE_KDTTREENODE_H

#include <optional>
#include <exception>
#include "Point2D.h"

namespace KDTree {

enum class NODESPLIT
{
    X, Y
};

/// \class KDTreeNode
/// \brief Simple tree leaf structure with alternating half splitting leafs ... for now.
/// \details Simple tree leaf structure
///          - Tree to incrementally add points to the structure.
///          - Get the closest point to a given input.
///          - Does not check for duplicates, expect unique points.
/// \tparam PointType the point type that the container will use.
/// \tparam NumVerticesInLeaf The number of points per leaf.
template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
class KDTreeNode
{
    using AccessX = Access<PointType, CoordinateType, 0>;
    using AccessY = Access<PointType, CoordinateType, 1>;
public:

    /// \brief Default constructor
    KDTreeNode() = default;

    /// \brief Custom constructor taking in the leaf nodes.
    /// \param one first leaf node.
    /// \param two second leaf node.
    KDTreeNode(std::unique_ptr<KDTreeNode>&& one, std::unique_ptr<KDTreeNode>&& two);

    /// \brief Custom constructor taking in the dimension of the leaf.
    /// \param min The minimum value of the leaf.
    /// \param extent The size of the leaf.
    /// \param split the direction the leaf is split for the children leafs.
    KDTreeNode(Point2D <CoordinateType> min, Point2D <CoordinateType> extent, NODESPLIT split) : min_(min),
                                                                                                 extent_(extent),
                                                                                                 split_(split) {}

    /// \brief Create a Parent leaf to contain the given point. This can create parents leafs recursively.
    /// \param point The point you want to insert into the domain that is out of the current domain.
    [[nodiscard]] std::unique_ptr<KDTreeNode> createParentLeaf(const PointType& point);

    /// \brief Check if the point is contained in the tree node.
    /// \param point The point to check if it is contained in the node.
    [[nodiscard]] bool isPointInCurrentDomain(const PointType& point);

    /// \brief Initialize the container to fit all the points.
    void initializeContainerDomainToFitAllPoints();

    /// \brief get the closest points to the input.
    /// \tparam T allows to get the closest point given any point representation having .x() and .y() functions.
    /// \param point The point you want to get the closest point from.
    template<typename T>
    void getClosestPoint(PointType& closestPoint, const T& point) const;

    /// \brief get the number of points in the container.
    [[nodiscard]] size_t size() { return numPoints_; }

    /// \brief Check if the container is empty.
    [[nodiscard]] bool empty() { return numPoints_ == 0; }

    /// \brief Insert a point in the storage.
    /// \param point The point to insert.
    template<typename T>
    void emplace(T&& point);

protected:
    std::array<std::unique_ptr<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>, 2> leafs_;
    std::array<PointType, NumVerticesInLeaf> points_;
    size_t numPoints_ = 0;
    std::optional<NODESPLIT> split_;
    std::optional<double> splitLocation_;
    std::optional<Point2D < CoordinateType>> min_;
    std::optional<Point2D < CoordinateType>> extent_;

    /// \brief Create children leafs and add the points into the children.
    void createLeafs();

    /// \brief Get the index for which leaf the point would be stored.
    [[nodiscard]] size_t getLeafIndex(const PointType& point) const;

    /// \brief Emplace the point in a child leaf.
    template<typename T>
    void emplaceInChildLeaf(T&& point);
};


template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::KDTreeNode(std::unique_ptr<KDTreeNode>&& one,
                                                                            std::unique_ptr<KDTreeNode>&& two)
{
    numPoints_ = one->numPoints_ + two->numPoints_;
    if (*one->split_ == NODESPLIT::X)
    {
        split_ = NODESPLIT::Y;
        extent_ = Point2D(one->extent_->x(), 2.0 * one->extent_->y());
    }
    else
    {
        split_ = NODESPLIT::X;
        extent_ = Point2D(2.0 * one->extent_->x(), one->extent_->y());
    }
    min_ = Point2D(std::min(one->min_->x(), two->min_->x()), std::min(one->min_->y(), two->min_->y()));
    if (*split_ == NODESPLIT::X)
    {
        if (one->min_->x() < two->min_->x())
        {
            leafs_[0] = std::move(one);
            leafs_[1] = std::move(two);
        }
        else
        {
            leafs_[0] = std::move(two);
            leafs_[1] = std::move(one);
        }
    }
    else
    {
        if (one->min_->y() < two->min_->y())
        {
            leafs_[0] = std::move(one);
            leafs_[1] = std::move(two);
        }
        else
        {
            leafs_[0] = std::move(two);
            leafs_[1] = std::move(one);
        }
    }
    if (*split_ == NODESPLIT::X)
    {
        splitLocation_ = min_->x() + 0.5 * extent_->x();
    }
    else
    {
        splitLocation_ = min_->y() + 0.5 * extent_->y();
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline std::unique_ptr<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>
KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::createParentLeaf(const PointType& point)
{
    if (*split_ == NODESPLIT::X)
    {
        Point2D min = *min_;
        if (AccessY::get(point) > min_->y())
        {
            min.y() += extent_->y();
        }
        else
        {
            min.y() -= extent_->y();
        }
        return std::make_unique<KDTreeNode>(min, *extent_, NODESPLIT::X);
    }
    else
    {
        Point2D min = *min_;
        if (AccessX::get(point) > min_->x())
        {
            min.x() += extent_->x();
        }
        else
        {
            min.x() -= extent_->x();
        }
        return std::make_unique<KDTreeNode>(min, *extent_, NODESPLIT::Y);
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline bool KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::isPointInCurrentDomain(const PointType& point)
{
    if (!min_)
    {
        return true;
    }
    return !(AccessX::get(point) < min_->x() || AccessY::get(point) < min_->y() ||
             AccessX::get(point) >= min_->x() + extent_->x() ||
             AccessY::get(point) >= min_->y() + extent_->y());

}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline void KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::initializeContainerDomainToFitAllPoints()
{
    Point2D<CoordinateType> min(AccessX::get(points_.front()), AccessY::get(points_.front()));
    Point2D<CoordinateType> max(AccessX::get(points_.front()), AccessY::get(points_.front()));
    Point2D<CoordinateType> extent(0, 0);
    for (auto& point : points_)
    {
        min.x() = std::min(min.x(), AccessX::get(point));
        min.y() = std::min(min.y(), AccessY::get(point));
        max.x() = std::max(max.x(), AccessX::get(point));
        max.y() = std::max(max.y(), AccessY::get(point));
    }
    extent.x() = max.x() - min.x();
    extent.y() = max.y() - min.y();
    extent.x() = std::max(extent.x() * 1.1, extent.y() * 1.1);
    extent.y() = extent.x();
    min_ = min;
    extent_ = extent;
    split_ = NODESPLIT::X;
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<typename T>
inline void KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::emplace(T&& point)
{
    if (!splitLocation_ && numPoints_ < points_.size())
    {
        points_[numPoints_] = std::forward<T>(point);
        numPoints_++;
    }
    else
    {
        if (!min_)
        {
            initializeContainerDomainToFitAllPoints();
        }
        if (!splitLocation_)
        {
            for (auto p : points_)
            {
                emplaceInChildLeaf(p);
            }
        }
        emplaceInChildLeaf(point);
        numPoints_++;
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<typename T>
inline void
KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::getClosestPoint(PointType& closestPoint, const T& point) const
{
    if (splitLocation_)
    {
        size_t leafIndex = getLeafIndex(point);
        leafs_[leafIndex]->getClosestPoint(closestPoint, point);
        double distSquared = (AccessX::get(point) - AccessX::get(closestPoint)) *
                             (AccessX::get(point) - AccessX::get(closestPoint)) +
                             (AccessY::get(point) - AccessY::get(closestPoint)) *
                             (AccessY::get(point) - AccessY::get(closestPoint));
        if ((*split_ == NODESPLIT::X &&
             distSquared >= (AccessX::get(point) - *splitLocation_) * (AccessX::get(point) - *splitLocation_)) ||
            (*split_ == NODESPLIT::Y &&
             distSquared >= (AccessY::get(point) - *splitLocation_) * (AccessY::get(point) - *splitLocation_)))
        {
            leafs_[1 - leafIndex]->getClosestPoint(closestPoint, point);
        }
    }
    else
    {
        double distSquared = (AccessX::get(point) - AccessX::get(closestPoint)) *
                             (AccessX::get(point) - AccessX::get(closestPoint)) +
                             (AccessY::get(point) - AccessY::get(closestPoint)) *
                             (AccessY::get(point) - AccessY::get(closestPoint));
        for (size_t i = 0; i < numPoints_; ++i)
        {
            const auto& p = points_[i];
            double distanceToPointSquared =
                    (AccessX::get(point) - AccessX::get(p)) * (AccessX::get(point) - AccessX::get(p)) +
                    (AccessY::get(point) - AccessY::get(p)) * (AccessY::get(point) - AccessY::get(p));
            if (distanceToPointSquared < distSquared)
            {
                distSquared = distanceToPointSquared;
                closestPoint = p;
            }
        }
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline void KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::createLeafs()
{
    double mid = 0;

    if (*split_ == NODESPLIT::X)
    {
        mid += min_->x() + 0.5 * extent_->x();
    }
    else
    {
        mid += min_->y() + 0.5 * extent_->y();
    }
    splitLocation_ = mid;
    if (*split_ == NODESPLIT::X)
    {
        leafs_[0] = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(*min_, Point2D((mid -
                                                                                                               min_->x()),
                                                                                                              extent_->y()),
                                                                                               NODESPLIT::Y);
        leafs_[1] = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(
                Point2D(min_->x() + (mid - min_->x()), min_->y()),
                Point2D(extent_->x() - (mid - min_->x()), extent_->y()), NODESPLIT::Y);
    }
    else
    {
        leafs_[0] = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(*min_,
                                                                                               Point2D(extent_->x(),
                                                                                                       (mid -
                                                                                                        min_->y())),
                                                                                               NODESPLIT::X);
        leafs_[1] = std::make_unique<KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>>(
                Point2D(min_->x(), min_->y() + (mid - min_->y())),
                Point2D(extent_->x(), extent_->y() - (mid - min_->y())), NODESPLIT::X);
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
inline size_t KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::getLeafIndex(const PointType& point) const
{
    if (*split_ == NODESPLIT::X)
    {
        return static_cast<size_t>(AccessX::get(point) > *splitLocation_);
    }
    else
    {
        return static_cast<size_t>(AccessY::get(point) > *splitLocation_);
    }
}

template<typename PointType, typename CoordinateType, size_t NumVerticesInLeaf>
template<typename T>
inline void KDTreeNode<PointType, CoordinateType, NumVerticesInLeaf>::emplaceInChildLeaf(T&& point)
{
    size_t i = getLeafIndex(point);
    if (!splitLocation_)
    {
        createLeafs();
    }
    leafs_.at(i)->emplace(std::forward<T>(point));
}

}

#endif // header guard