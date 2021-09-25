/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at https://mozilla.org/MPL/2.0/.
/// author: Andre Fecteau <andre.fecteau1@gmail.com>

#define CATCH_CONFIG_MAIN

#include "../lib/catch/include/catch.hpp"
#include "../include/KDTree.h"


TEST_CASE("emplace(x,y)", "test_KDTree")
{
    KDTree::KDTree tree;
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = 0; j < 10; ++j)
        {
            tree.emplace(i, j);
        }
    }
    REQUIRE(tree.size() == 100);
}

TEST_CASE("emplace<T>(T)", "test_KDTree")
{
    KDTree::KDTree tree;
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = 0; j < 10; ++j)
        {
            tree.emplace(KDTree::Point2D(static_cast<double>(i), static_cast<double>(j)));
        }
    }
    REQUIRE(tree.size() == 100);
}


TEST_CASE("Small tree test", "closest point")
{
    srand(time(nullptr));
    KDTree::KDTree tree;
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = 0; j < 10; ++j)
        {
            bool inserted = false;
            while (!inserted)
            {
                size_t x = (rand() % 10000);
                size_t y = (rand() % 10000);
                if (!tree.empty())
                {
                    auto closestPoint = tree.getClosestPoint(KDTree::Point2D(x * 1e-5, y * 1e-5));
                    if (std::abs(closestPoint.x() - x * 1e-5) > 1e-5 ||
                        std::abs(closestPoint.y() - y * 1e-5) > 1e-5)
                    {
                        tree.emplace(x * 1e-5, y * 1e-5);
                        inserted = true;
                    }
                }
                else
                {
                    tree.emplace(x * 1e-5, y * 1e-5);
                    inserted = true;
                }
            }
        }
    }
    REQUIRE(tree.size() == 100);
}

struct OtherPoint2D
{
    /// Default Constructor
    OtherPoint2D() = default;

    OtherPoint2D(const double& x, const double& y) : x(x), y(y) {}

    double x; ///< first dimension value
    double y; ///< second dimension value
};

namespace KDTree {
template<>
struct Access<OtherPoint2D, double, 0>
{
    static void set(OtherPoint2D& point, const double& value) { point.x = value; }

    static const double& get(const OtherPoint2D& point) { return point.x; }
};

template<>
struct Access<OtherPoint2D, double, 1>
{
    static void set(OtherPoint2D& point, const double& value) { point.y = value; }

    static const double& get(const OtherPoint2D& point) { return point.y; }
};
}

TEST_CASE("TemplateParameters", "PointType")
{
    KDTree::KDTree<OtherPoint2D> tree;
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = 0; j < 10; ++j)
        {
            tree.emplace(static_cast<double>(i), static_cast<double>(j));
        }
    }
    REQUIRE(tree.size() == 100);
}