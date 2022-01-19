#include "color.h"


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



Color& Color::from_CIE_RGB(const Vec& rgb) {
    x = 0.4124 * rgb.x + 0.3576 * rgb.y + 0.1805 * rgb.z;
    y = 0.2126 * rgb.x + 0.7152 * rgb.y + 0.0722 * rgb.z;
    z = 0.0193 * rgb.x + 0.1192 * rgb.y + 0.9505 * rgb.z;
    return *this;
}


Color& Color::from_sRGB(Vec srgb) {
    srgb.x = gamma_expand(srgb.x);
    srgb.y = gamma_expand(srgb.y);
    srgb.z = gamma_expand(srgb.z);
    return from_CIE_RGB(srgb);
}


/* inline */ Color& Color::operator+=(const Color& vb) { 
    x += vb.x;
    y += vb.y;
    z += vb.z;
    return *this;
}


/* inline */ Color Color::operator+(const Color& vb) const { 
    Color res(*this);
    res += vb;
    return res;
}


/* inline */ Color& Color::operator-=(const Color& vb) {
    x -= vb.x;
    y -= vb.y;
    z -= vb.z;
    return *this;
}


/* inline */ Color Color::operator-(const Color& vb) const { 
    Color res(*this);
    res -= vb;
    return res;
}


/* inline */ Color& Color::operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
}


/* inline */ Color Color::operator/(double scale) const {
    Color res(*this);
    res /= scale;
    return res;
}


/* inline */ Color& Color::operator*=(double scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;   
}



/* inline */ Color operator*(Color lhs, double scale) {
    lhs *= scale;
    return lhs;
}


/* inline */ Color operator*(double scale, Color rhs) {
    rhs *= scale;
    return rhs;
}