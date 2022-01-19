#include "vector.h"

#include <cmath>
#include <stdio.h>
#include <stdexcept>

/* inline */ Vec& Vec::operator+=(const Vec& vb) { 
    x += vb.x;
    y += vb.y;
    z += vb.z;
    return *this;
}


/* inline */ Vec Vec::operator+(const Vec& vb) const { 
    Vec res(*this);
    res += vb;
    return res;
}


/* inline */ Vec& Vec::operator-=(const Vec& vb) {
    x -= vb.x;
    y -= vb.y;
    z -= vb.z;
    return *this;
}


/* inline */ Vec Vec::operator-(const Vec& vb) const { 
    Vec res(*this);
    res -= vb;
    return res;
}


/* inline */ Vec& Vec::operator%=(const Vec& vb) {
    double new_x = y * vb.z - z * vb.y;
    double new_y = z * vb.x - x * vb.z;
    double new_z = x * vb.y - y * vb.x;
    x = new_x;
    y = new_y;
    z = new_z;
    return *this;
}


/* inline */ Vec Vec::operator%(const Vec& vb) const { 
    Vec res(*this);
    res %= vb;
    return res;
}


/* inline */ Vec& Vec::operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
}


/* inline */ Vec Vec::operator/(double scale) const {
    Vec res(*this);
    res /= scale;
    return res;
}


/* inline */ Vec& Vec::operator*=(double scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;   
}


/* inline */ double Vec::norm() const {
    return std::sqrt(x*x + y*y + z*z);
}


/* inline */ Vec& Vec::normalize() {
    return *this *= 1 / this->norm();
}


/* inline */ Vec operator*(Vec lhs, double scale) {
    lhs *= scale;
    return lhs;
}


/* inline */ Vec operator*(double scale, Vec rhs) {
    rhs *= scale;
    return rhs;
}


/* inline */ double& Vec::operator[](int idx) {
    switch (idx) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("0: x, 1: y, 2: z");
    }
}

/* inline */ double Vec::operator[](int idx) const {
    switch (idx) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("0: x, 1: y, 2: z");
    }
}


/* inline */ double Vec::dot(const Vec& vb) const {
    return x * vb.x + y * vb.y + z * vb.z;
}


void Vec::show() const {
    printf("[*] x = %f, y = %f, z = %f\n", x, y, z);
    printf("[*] ||v|| = %f\n", this->norm());
}