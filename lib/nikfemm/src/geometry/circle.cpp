#include <math.h>
#include <assert.h>
#include <algorithm>
#include <random>
#include <vector>
#include <unordered_set>

#include <constants.hpp>

#include "circle.hpp"

namespace nikfemm {
    Circle::Circle(Vector center, double radius) {
        this->center = center;
        this->radius = radius;
    }

    Circle::Circle() {
        this->center = Vector();
        this->radius = 0;
    }

    bool Circle::operator==(const Circle& c) const {
        return (center == c.center) && (radius == c.radius);
    }

    bool Circle::operator!=(const Circle& c) const {
        return !(*this == c);
    }
    
    Circle Circle::trivialCircleFromPoints(std::vector<Vector> points) {
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

    Circle Circle::welzlHelper(std::vector<Vector> points, std::vector<Vector> R, uint32_t n) {
        if (n == 0 || R.size() == 3) {
            return trivialCircleFromPoints(R);
        }

        uint32_t k = rand() % n;
        Vector p = points[k];

        std::swap(points[k], points[n - 1]);
        Circle c = welzlHelper(points, R, n - 1);

        if (c.contains(p)) {
            return c;
        }

        R.push_back(p);
        return welzlHelper(points, R, n - 1);
    }

    Circle Circle::getMinimumEnclosingCircle(std::vector<Vector> points) {
        auto rng = std::default_random_engine {};
        std::shuffle(points.begin(), points.end(), rng);
        return welzlHelper(points, std::vector<Vector>(), points.size());
    }

    Circle Circle::getMinimumEnclosingCircle(std::unordered_set<Vector> points) {
        std::vector<Vector> p;
        p.reserve(points.size());
        for (auto it = points.begin(); it != points.end(); ) {
            p.push_back(std::move(points.extract(it++).value()));
        }
        return getMinimumEnclosingCircle(p);
    }

    Circle Circle::getEnclosingCircle(std::vector<Vector> points) {
        double maxmod = 0;
        Vector center = {0, 0};
        for (auto p : points) {
            double mod = p.norm();
            if (mod > maxmod) {
                maxmod = mod;
            }
            center += p;
        }
        center /= points.size();
        return Circle(center, maxmod);
    }
}