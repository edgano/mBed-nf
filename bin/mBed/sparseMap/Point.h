#ifndef POINT_H
#define POINT_H

#include<vector>
#include <iostream>
#include "types.h"
#ifdef DEBUG_MODE
#define ASSERT(condition) if (!(condition)) assertionFailure(__FILE__, __LINE__);
#else
#define ASSERT(condition)
#endif

using namespace std;

/**
 * Point class
 **/
class Point {
public:
    Point() {}
    Point(int numDimensions) {myData = vector<Scalar>(numDimensions, 0.0); myDim = numDimensions /*FS*/;}
    Point(const Point &rhs) {myData = rhs.myData;}

    inline int getNumDimensions() const {return (int)myData.size();}
    inline Scalar operator[](int index) const {ASSERT(index >= 0 && index < getNumDimensions()); return myData[index];}
    inline Scalar &operator[](int index) {ASSERT(index >= 0 && index < getNumDimensions()); return myData[index];}

    int myDim/*FS*/;

private:
	vector<Scalar> myData;
};


#endif
