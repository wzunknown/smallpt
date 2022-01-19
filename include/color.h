#pragma once
#include "vector.h"
#include <cmath>


using Color_8bit = unsigned char;

/**
 * @brief CIE XYZ color space
 * 
 */
class Color: public Vec {
public:
    Color(double x_=0, double y_=0, double z_=0) : Vec(x_, y_, z_) {}

    Vec to_CIE_RGB();
    Vec to_sRGB();
    Color& from_CIE_RGB(const Vec& rgb);
    Color& from_sRGB(Vec srgb);


    Color& operator+=(const Color& vb);
    Color operator+(const Color& vb) const;
    Color& operator-=(const Color& vb);
    Color operator-(const Color& vb) const;
    Color& operator*=(double scale);
    Color& operator/=(double scale);
    Color operator/(double scale) const;


    Color mult(const Color& color) const {
        return Color(x * color.x, y * color.y, z * color.z);
    }


    static inline double clamp(double val) {
        return val < 0? 0 : val > 1? 1 : val;
    }


    static inline Color_8bit to_8bit(double val) {
        // int res = (int) round(clamp(val) * 255);
        // return res < 0? 0 : res > 255? 255 : (Color_8bit)res;
        return (Color_8bit) round(clamp(val) * 255);
    }

    static inline double from_8bit(Color_8bit val) {
        return val / 255.;
    } 
};


extern Color operator*(Color lhs, double scale);
extern Color operator*(double scale, Color rhs);