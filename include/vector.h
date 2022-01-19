#pragma once


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


    double& operator[](int idx);
    double operator[](int idx) const;


    
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


extern Vec operator*(Vec lhs, double scale);
extern Vec operator*(double scale, Vec rhs);
