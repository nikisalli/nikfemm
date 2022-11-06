#include <math.h>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <unordered_set>

#include <constants.hpp>

#include "circle.hpp"

namespace nikfemm {
    Circle::Circle(Point center, double radius) {
        this->center = center;
        this->radius = radius;
    }

    Circle::Circle() {
        this->center = Point();
        this->radius = 0;
    }

    bool Circle::operator==(const Circle& c) const {
        return fabs(this->radius - c.radius) < EPSILON && fabs(this->center.x - c.center.x) < EPSILON && fabs(this->center.y - c.center.y) < EPSILON;
    }

    bool Circle::operator!=(const Circle& c) const {
        return !(*this == c);
    }
    
    Circle Circle::trivialCircleFromPoints(std::vector<Point> points) {
        assert(points.size() <= 3);
        if (points.size() == 0) {
            return Circle();
        } else if (points.size() == 1) {
            return Circle(points[0], 0);
        } else if (points.size() == 2) {
            return getCircleFromPoints(points[0], points[1]);
        }

        // check if minimum enclosing circle can be found with 2 points
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
    
                Circle c = getCircleFromPoints(points[i], points[j]);
                if (c.containsPoints(points)) {
                    return c;
                }
            }
        }

        return getCircleFromPoints(points[0], points[1], points[2]);
    }

    Circle Circle::welzlHelper(std::vector<Point> points, std::vector<Point> R, uint32_t n) {
        if (n == 0 || R.size() == 3) {
            return trivialCircleFromPoints(R);
        }

        uint32_t k = rand() % n;
        Point p = points[k];

        std::swap(points[k], points[n - 1]);
        Circle c = welzlHelper(points, R, n - 1);

        if (c.contains(p)) {
            return c;
        }

        R.push_back(p);
        return welzlHelper(points, R, n - 1);
    }

    Circle Circle::getMinimumEnclosingCircle(std::vector<Point> points) {
        std::random_shuffle(points.begin(), points.end());
        return welzlHelper(points, std::vector<Point>(), points.size());
    }

    Circle Circle::getMinimumEnclosingCircle(std::unordered_set<Point> points) {
        std::vector<Point> p;
        p.reserve(points.size());
        for (auto it = points.begin(); it != points.end(); ) {
            p.push_back(std::move(points.extract(it++).value()));
        }
        return getMinimumEnclosingCircle(p);
    }
}