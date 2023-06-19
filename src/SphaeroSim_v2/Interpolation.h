#ifndef _SPHAEROSIM_INTERPOLATION_H_
#define _SPHAEROSIM_INTERPOLATION_H_

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_
#include "EigenLibrary.h"
#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

namespace SphaeroSim {

    enum InterpolationType {
        kLinear,
        kSpline
    };

    template<class T>
    class Interpolation {
    public:
        explicit Interpolation(const InterpolationType type) {
            type_ = type;
        }

        typedef std::pair<const double, const T> ref_point;
        typedef std::pair<const double, const T *> ref_point_ptr;

        inline
        T ExecuteInterpolation(const ref_point &point_1,
                               const ref_point &point_2,
                               const ref_point_ptr &before,
                               const ref_point_ptr &after,
                               const double delta) const {
            if (GetType() == kLinear) {
                return ExecuteLinearInterpolation(point_1, point_2, delta);
            } else if (GetType() == kSpline) {
                T ctrl_points[4];
                ControlPoints(point_1,
                              point_2,
                              before,
                              after,
                              ctrl_points);
                double relative_delta = (delta - point_1.first) / (point_2.first - point_1.first);
                return ExecuteSplineInterpolation(ctrl_points, relative_delta);
            } else {
                throw Exception("ExecuteInterpolation",
                                "Unsupported interpolation type.");
            }
        }

        // Getter
        inline const InterpolationType GetType() const {
            return type_;
        }

    private:
        T ExecuteLinearInterpolation(const ref_point &point_1,
                                     const ref_point &point_2,
                                     const double delta) const {
            double interval = point_2.first - point_1.first;
            T slope = (point_2.second - point_1.second) / interval;
            T result = point_1.second + slope * (delta - point_1.first);
            return result;
        }

        T ExecuteSplineInterpolation(const T ctrl_pts[4],
                                     const double relative_delta) const {
            double rInv = 1.0 - relative_delta;
            T res = rInv * rInv * rInv * ctrl_pts[0] +
                    3.0 * relative_delta * rInv * rInv * ctrl_pts[1] +
                    3.0 * relative_delta * relative_delta * rInv * ctrl_pts[2] +
                    relative_delta * relative_delta * relative_delta * ctrl_pts[3];
            return res;
        }

        const T slope(const ref_point &point,
                      const ref_point &prev,
                      const ref_point &next) const {
            T slope1 = (point.second - prev.second) / (point.first - prev.first);
            T slope2 = (next.second - point.second) / (next.first - point.first);

            return 0.5 * (slope1 + slope2);
        }

        void ControlPoints(const ref_point &point_1,
                           const ref_point &point_2,
                           const ref_point_ptr &before,
                           const ref_point_ptr &after,
                           T dst_points[4]) const {
            T slope_pos;
            T slope_neg;

            dst_points[0] = point_1.second;
            dst_points[3] = point_2.second;

            if (before.second == NULL) {
                slope_pos = (point_2.second - point_1.second) /
                            (point_2.first - point_1.first);
            } else {
                slope_pos = slope(point_1,
                                  ref_point(point_1.first - before.first,
                                            *before.second),
                                  point_2);
            }

            if (after.second == NULL) {
                slope_neg = (point_2.second - point_1.second) /
                            (point_2.first - point_1.first);
            } else {
                slope_neg = slope(point_2,
                                  point_1,
                                  ref_point(point_2.first + after.first,
                                            *after.second));
            }

            dst_points[1] = slope_pos * (point_2.first - point_1.first) / 3.0 +
                            point_1.second;
            dst_points[2] = -slope_neg * (point_2.first - point_1.first) / 3.0 +
                            point_2.second;
        }

        InterpolationType type_;
    };
}  // namespace SphaeroSim
#endif  // _SPHAEROSIM_INTERPOLATION_H_