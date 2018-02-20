#include "PointUtils.h"

/**
 * Outputs a point in human readable form.
 **/
ostream &operator<<(ostream &out, const Point &rhs) 
{
    out << "(";
    for (int i = 0; i < rhs.getNumDimensions(); i++)
        out << rhs[i] << (i+1 == rhs.getNumDimensions()? ")" : ", ");
    return out;
}

/**
 * Outputs a list of points in human readable form.
 **/
ostream &operator<<(ostream &out, const vector<Point> &rhs) 
{
    out << "{";
    for (int i = 0; i < (int)rhs.size(); i++)
        out << rhs[i] << (i+1 == int(rhs.size())? "}" : ",\n ");
    return out;
}

/**
 * Returns the distance squared between two points.
 **/
/*Scalar distSq(const Point &lhs, const Point &rhs) 
{
    ASSERT(lhs.getNumDimensions() == rhs.getNumDimensions());
    Scalar result = 0;
    for (int i = 0; i < lhs.getNumDimensions(); i++)
        result += (lhs[i] - rhs[i]) * (lhs[i] - rhs[i]);
    return result;
}*/

/**
 * Reflects this point about the origin and returns the result.
 **/
Point operator-(const Point &x) 
{
    Point p(x.getNumDimensions());
    for (int i = 0; i < x.getNumDimensions(); i++)
        p[i] = -x[i];
    return p;
}

/**
 * Adds two points as vectors.
 **/
Point operator+(const Point &lhs, const Point &rhs) 
{
    Point p = lhs;
    ASSERT(lhs.getNumDimensions() == rhs.getNumDimensions());
    
    for (int i = 0; i < lhs.getNumDimensions(); i++)
        p[i] += rhs[i];
    return p;
}

/**
 * Subtracts two points as vectors.
 **/
Point operator-(const Point &lhs, const Point &rhs) 
{
    return lhs + (-rhs);
}

/**
 * Scales this point as a vector.
 **/
Point operator*(const Point &lhs, Scalar f) 
{
    Point p = lhs;
    for (int i = 0; i < lhs.getNumDimensions(); i++)
        p[i] *= f;
    return p;
}

Point operator*(Scalar f, const Point &lhs) 
{
    return lhs * f;
}

/**
 * Scales this point as a vector.
 **/
Point operator/(const Point &lhs, Scalar f) 
{
    ASSERT(f != 0);
    return lhs * (1/f);
}
