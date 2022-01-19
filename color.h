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
};

static double gamma_compress(double linear) {
    return linear <= 0.0031308? linear * 12.92 : 1.055 * pow(linear, 1.0 / 2.4) - 0.055;
}

static double gamma_expand(double non_linear) {
    return non_linear <= 0.04045? non_linear / 12.92 : pow((non_linear + 0.055)/1.055, 2.4); 
}


Vec Color::to_CIE_RGB() {
    return Vec(
        3.2406 * x - 1.5372 * y - 0.4986 * z,
        -0.9689 * x + 1.8758 * y + 0.0415 * z,
        0.0557 * x - 0.2040 * y + 1.0570 * z
    );
    // return Vec(x, y, z);
}


Vec Color::to_sRGB() {
    Vec srgb = to_CIE_RGB();
    srgb.x = clamp(gamma_compress(srgb.x));
    srgb.y = clamp(gamma_compress(srgb.y));
    srgb.z = clamp(gamma_compress(srgb.z));
    return srgb;
}


inline Color& Color::operator+=(const Color& vb) { 
    x += vb.x;
    y += vb.y;
    z += vb.z;
    return *this;
}


inline Color Color::operator+(const Color& vb) const { 
    Color res(*this);
    res += vb;
    return res;
}


inline Color& Color::operator-=(const Color& vb) {
    x -= vb.x;
    y -= vb.y;
    z -= vb.z;
    return *this;
}


inline Color Color::operator-(const Color& vb) const { 
    Color res(*this);
    res -= vb;
    return res;
}


inline Color& Color::operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
}


inline Color Color::operator/(double scale) const {
    Color res(*this);
    res /= scale;
    return res;
}


inline Color& Color::operator*=(double scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;   
}



inline Color operator*(Color lhs, double scale) {
    lhs *= scale;
    return lhs;
}


inline Color operator*(double scale, Color rhs) {
    rhs *= scale;
    return rhs;
}