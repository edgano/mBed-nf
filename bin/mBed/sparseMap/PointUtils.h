#ifndef POINTUTILS_H
#define POINTUTILS_H

#include<vector>
#include <iostream>
#include "Point.h"
#include "types.h"
/**
 * Point utilities
 **/
ostream &operator<<(ostream &out, const Point &rhs);
ostream &operator<<(ostream &out, const vector<Point> &rhs);
//inline Scalar distSq(const Point &lhs, const Point &rhs);
Point operator-(const Point &x);
Point operator+(const Point &lhs, const Point &rhs);
Point operator-(const Point &lhs, const Point &rhs);
Point operator*(const Point &lhs, Scalar f);
Point operator*(Scalar f, const Point &lhs);
Point operator/(const Point &lhs, Scalar f);

#endif
