#pragma once
#include <cmath>
#include <stdio.h>


// vector operations
class Vec {
public:
    double x, y, z;
    Vec(double x_=0, double y_=0, double z_=0) : x{x_}, y{y_}, z{z_} {}

    Vec& operator+=(const Vec& vb);
    Vec operator+(const Vec& vb) const;
    
    Vec& operator-=(const Vec& vb);
    Vec operator-(const Vec& vb) const;


    /**
     * @brief cross product of vectors
     * 
     * @param vb 
     * @return Vec& 
     */
    Vec& operator%=(const Vec& vb);
    Vec operator%(const Vec& vb) const;

    Vec& operator/=(double scale);
    Vec operator/(double scale) const;


    /**
     * @brief scalar product
     * 
     * @param scale 
     * @return Vec& 
     */
    Vec& operator*=(double scale);


    
    /**
     * @brief dot product of vectors
     * 
     * @param vb 
     * @return double 
     */
    double dot(const Vec& vb) const;

    double norm() const;
    Vec& normalize();
    void show() const;
};


inline Vec& Vec::operator+=(const Vec& vb) { 
    x += vb.x;
    y += vb.y;
    z += vb.z;
    return *this;
}


inline Vec Vec::operator+(const Vec& vb) const { 
    Vec res(*this);
    res += vb;
    return res;
}


inline Vec& Vec::operator-=(const Vec& vb) {
    x -= vb.x;
    y -= vb.y;
    z -= vb.z;
    return *this;
}


inline Vec Vec::operator-(const Vec& vb) const { 
    Vec res(*this);
    res -= vb;
    return res;
}


inline Vec& Vec::operator%=(const Vec& vb) {
    double new_x = y * vb.z - z * vb.y;
    double new_y = z * vb.x - x * vb.z;
    double new_z = x * vb.y - y * vb.x;
    x = new_x;
    y = new_y;
    z = new_z;
    return *this;
}


inline Vec Vec::operator%(const Vec& vb) const { 
    Vec res(*this);
    res %= vb;
    return res;
}


inline Vec& Vec::operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
}


inline Vec Vec::operator/(double scale) const {
    Vec res(*this);
    res /= scale;
    return res;
}


inline Vec& Vec::operator*=(double scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;   
}


inline double Vec::norm() const {
    return sqrt(x*x + y*y + z*z);
}


inline Vec& Vec::normalize() {
    return *this *= 1 / this->norm();
}


inline Vec operator*(Vec lhs, double scale) {
    lhs *= scale;
    return lhs;
}


inline Vec operator*(double scale, Vec rhs) {
    rhs *= scale;
    return rhs;
}


inline double Vec::dot(const Vec& vb) const {
    return x * vb.x + y * vb.y + z * vb.z;
}


void Vec::show() const {
    printf("[*] x = %f, y = %f, z = %f\n", x, y, z);
    printf("[*] ||v|| = %f\n", this->norm());
}