#ifndef UTILS_SPARSEMAP_H
#define UTILS_SPARSEMAP_H


/**
 * Includes
 **/
#include <iostream>
#include <vector>
#include <memory>
#include "Point.h"
#include "PointUtils.h"
#include "types.h"
#include "../tree/UPGMA/Node.h"
using namespace std;

/**
 * Type definitions
 **/

struct KMeansResult
{
    auto_ptr<vector<Point> > centers;
    auto_ptr<vector<int> > clusterAssignment;
};

/**
 * General utilities
 **/
//#define DEBUG_MODE
#ifdef DEBUG_MODE
#define ASSERT(condition) if (!(condition)) assertionFailure(__FILE__, __LINE__);
#else
#define ASSERT(condition)
#endif
void fatalError(const char *message);
void assertionFailure(const char *filename, int lineNumber);

/**
 * Randomization utilities
 **/
Scalar getRandomScalar();
Point getRandomPoint(int numDimensions);
Point getByNormDist(const Point &center, Scalar variance);

/**
 * K-means utilities
 **/
//vector<Point> chooseUniformCenters(const vector<Point> &data, int numCenters);
auto_ptr<vector<Point> > chooseSmartCenters(const vector<Point>* data, int numCenters, int numLocalTries);
auto_ptr<KMeansResult> runKMeans(const vector<Point>* data, vector<Point> *initCenters, int &numIterations);
void printKMeansResults(KMeansResult* results);
Scalar euclideanDistance(Point* si, Point* sj);
//Scalar getKMeansPotential(const vector<Point> &data, const vector<Point> &centers);

void preorder(clustalw::Node* root);

#endif

