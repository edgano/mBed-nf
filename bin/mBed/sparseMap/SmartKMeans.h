#ifndef SMARTKMEANS_H
#define SMARTKMEANS_H

#include <vector>
#include <memory>
#include "Point.h"
#include "Utils.h"
/**
 * K-means with smart clustering, no retries at each stage.
 **/
class SmartKMeans
{
    public:
        string getName() const 
        {
            return "K-MEANS++";
        }
        auto_ptr<vector<Point> > initClustering(const vector<Point>* data, int numCenters) const 
        {
            return chooseSmartCenters(data, numCenters, 1);
        }
	auto_ptr<KMeansResult> doClustering(const vector<Point>* data, 
	                                        vector<Point>* clustering, int &numIterations) const 
        {
            return runKMeans(data, clustering, numIterations);
        }
};

#endif


