/**
 * Includes
 **/
#include "Utils.h"
#include "PointUtils.h"
#include <cmath>
#include <cstdlib>

inline Scalar distSq(const Point &lhs, const Point &rhs);

Scalar distSq(const Point &lhs, const Point &rhs) 
{
    ASSERT(lhs.getNumDimensions() == rhs.getNumDimensions());
    Scalar result = 0;
    for (int i = 0; i < lhs.getNumDimensions(); i++)
        result += (lhs[i] - rhs[i]) * (lhs[i] - rhs[i]);
    return result;
}


/**
 * Generates a fatal error resulting from an assertion failing at the given location.
 **/
void assertionFailure(const char *filename, int lineNumber) 
{
    //const char *errorMessage = "Assertion failed on line %d of file %s.";
    //char *buffer = new char[_scprintf(errorMessage, lineNumber, filename)+1];
    //sprintf(buffer, errorMessage, lineNumber, filename);
    //fatalError(buffer);
}




/////////////////////////////
// RANDOMIZATION UTILITIES //
/////////////////////////////

/**
 * Returns a number uniformly at random from [0,1).
 **/
Scalar getRandomScalar() 
{
    return Scalar(rand()) / RAND_MAX;
}

/**
 * Returns a point with coordinates chosen uniformly at random from [0,1).
 **/
Point getRandomPoint(int numDimensions) 
{
    Point p(numDimensions);
    for (int i = 0; i < numDimensions; i++)
        p[i] = getRandomScalar();
    return p;
}

/**
 * Returns a point chosen according to a given normal distribution.
 **/
Scalar normRange = 50;
int normIntervals = 50000;
vector<Scalar> stdNormCumDist;

void initStdNormCumDist() 
{
    stdNormCumDist.push_back(0);
    for (int i = 0; i <= normIntervals; i++) 
    {
        Scalar x = -normRange/2 + (normRange*i) / (normIntervals-1);
        Scalar y = exp(-(x*x)/2);
        stdNormCumDist.push_back(stdNormCumDist[i] + y);
    }
}

Scalar getByStdNormDist() 
{
    if (stdNormCumDist.size() == 0)
        initStdNormCumDist();
    
    Scalar maxY = stdNormCumDist[normIntervals];
    Scalar y = getRandomScalar() * maxY;

    // Binary search so that lb = ub
    int lb = 0, ub = normIntervals-1;
    while (lb < ub) 
    {
        int mb = (lb+ub+1)/2;
        if (stdNormCumDist[mb] > y) ub = mb-1;
        else lb = mb;
    }

    // Find an exact value by interpolation - it IS possible for division by 0 due to precision problems
    // so avoid that with a hack
    Scalar x1 = -normRange/2 + (normRange*lb) / (normIntervals-1);
    Scalar x2 = -normRange/2 + (normRange*(lb+1)) / (normIntervals-1);
    if (stdNormCumDist[lb] - stdNormCumDist[lb+1] == 0) return x1;
    return x1 + (x2-x1) * (y - stdNormCumDist[lb]) / (stdNormCumDist[lb+1] - stdNormCumDist[lb]);
}

Point getByNormDist(const Point &center, Scalar variance) 
{
    Point result = center;
    for (int i = 0; i < center.getNumDimensions(); i++)
        result[i] += getByStdNormDist() * variance;
    return result;
}

///////////////////////
// K-MEANS UTILITIES //
///////////////////////

/**
 * Given fixed data points and fixed centers, this returns the index of which center is closest
 * to each point.
 **/
auto_ptr<vector<int> > assignPointsToClusters(const vector<Point>* data, const vector<Point>* centers) 
{
    auto_ptr<vector<int> > result;
    result.reset(new vector<int>(data->size()));
    
    for (int i = 0; i < (int)data->size(); i++) 
    {
        Scalar bestValue = -1;
        int bestIndex;
        for (int j = 0; j < (int)centers->size(); j++) 
        {
            Scalar thisValue = distSq((*data)[i], (*centers)[j]);
            if (bestValue < 0 || thisValue < bestValue) 
            {
                bestValue = thisValue;
                bestIndex = j;
            }
        }
        (*result)[i] = bestIndex;
    }
    return result;
}

/**
 * Given fixed data points and a fixed cluster assignment, this returns the optimal centers for this
 * cluster assignment.
 **/
auto_ptr<vector<Point> > computeCenters(const vector<Point>* data, const vector<int> *clusterAssignment, int numClusters) 
{
    int i;
    vector<int> clusterSize(numClusters);
    vector<Point> clusterPosition(numClusters, Point((*data)[0].getNumDimensions()));

    for (i = 0; i < (int)clusterAssignment->size(); i++) 
    {
        clusterSize[(*clusterAssignment)[i]]++;    
        clusterPosition[(*clusterAssignment)[i]] = clusterPosition[(*clusterAssignment)[i]] + (*data)[i];
    }

    auto_ptr<vector<Point> > result;
    result.reset(new vector<Point>);
    for (i = 0; i < numClusters; i++)
        if (clusterSize[i] > 0)
            result->push_back(clusterPosition[i] / clusterSize[i]);

    return result;
}

/**
 * Chooses a number of centers uniformly at random from the data set.
 **/
/*vector<Point> chooseUniformCenters(const vector<Point> &data, int numCenters) 
{
    ASSERT( numCenters > 0 && numCenters <= (int)data.size() );
    int i;

    // Create a list of unchosen centers
    vector<Point> centers;
    vector<int> centerIndices;
    for (i = 0; i < (int)data.size(); i++)
        centerIndices.push_back(i);

    // Choose each center one at a time, keeping the list up to date
    for (i = 0; i < numCenters; i++) 
    {
        int index = (int)(getRandomScalar() * centerIndices.size());
        centers.push_back( data[centerIndices[index]] );
        centerIndices[index] = centerIndices[ int(centerIndices.size()) - 1 ];
        centerIndices.pop_back();
    }

    return centers;
}*/

/**
 * Chooses a number of centers from the data set as follows:
 *  - One center is chosen randomly.
 *  - Now repeat numCenters-1 times:
 *      - Repeat numLocalTries times:
 *          - Add a point x with probability proportional to the distance squared from x
 *            to the closest existing center
 *      - Add the point chosen above that results in the smallest potential.
 **/
auto_ptr<vector<Point> > chooseSmartCenters(const vector<Point>* data, int numCenters, int numLocalTries) 
{
    ASSERT( numCenters > 0 && numCenters <= (int)data->size() );
    int i;
    Scalar currentPot = 0;
    auto_ptr<vector<Point> > centers;
    centers.reset(new vector<Point>);
    vector<Scalar> closestDistSq;

    // Choose one random center and set the closestDistSq values
    int index = (int)(getRandomScalar() * data->size());
    centers->push_back((*data)[index]);
    for (i = 0; i < (int)data->size(); i++) 
    {
        closestDistSq.push_back( distSq((*data)[i], (*data)[index]) );
        currentPot += closestDistSq[i];
    }

    // Choose each center
    for (int centerCount = 1; centerCount < numCenters; centerCount++) 
    {

        // Repeat several trials
        Scalar bestNewPot = -1;
        int bestNewIndex;
        for (int localTrial = 0; localTrial < numLocalTries; localTrial++) 
        {
    
        // Choose our center - have to be slightly careful to return a valid 
        // answer even accounting
        // for possible rounding errors
            Scalar randVal = getRandomScalar() * currentPot;
            for (index = 0; index < (int)data->size()-1; index++) 
            {
                if (randVal <= closestDistSq[index])
                    break;
                else
                    randVal -= closestDistSq[index];
            }

            // Compute the new potential
            Scalar newPot = 0;
            for (i = 0; i < (int)data->size(); i++)
                newPot += min( distSq((*data)[i], (*data)[index]), closestDistSq[i] );

            // Store the best result
            if (bestNewPot < 0 || newPot < bestNewPot) 
            {
                bestNewPot = newPot;
                bestNewIndex = index;
            }
        }

        // Add the appropriate center
        centers->push_back((*data)[bestNewIndex]);
        currentPot = bestNewPot;
        for (i = 0; i < (int)data->size(); i++)
            closestDistSq[i] = min( distSq((*data)[i], (*data)[bestNewIndex]),
                                           closestDistSq[i] );
    }

    return centers;
}

/**
 * Runs k-means on the given data points with the given initial centers.
 * It iterates up to numIterations times, and for each iteration it does, it decreases numIterations by 1.
 **/
auto_ptr<KMeansResult> runKMeans(const vector<Point>* data, vector<Point>* initCenters, int &numIterations) 
{
    cout << "Num of iterations " << numIterations << "\n";
    auto_ptr<vector<Point> > centers;
    centers.reset(new vector<Point>(*initCenters));
    initCenters->clear();
    
    auto_ptr<vector<int> > clusterAssignment = assignPointsToClusters(data, centers.get());

    while (numIterations > 0) 
    {
        centers = computeCenters(data, clusterAssignment.get(), (int)centers->size());
        auto_ptr<vector<int> > newClusterAssignment = 
                                            assignPointsToClusters(data, centers.get());
        bool isProgress = false;
        for (int i = 0; i < (int)clusterAssignment->size(); i++)
            if ((*clusterAssignment)[i] != (*newClusterAssignment)[i])
                isProgress = true;

        if (!isProgress)
            break;

        clusterAssignment = newClusterAssignment;
        numIterations--;
    }
    
    //cout << "Here are the cluster members\n";
    //for(int i = 0; i < clusterAssignment.size(); i++)
    //{
    //    cout << clusterAssignment[i] << " ";
    //}
    //cout << "\n";
    auto_ptr<KMeansResult> results;
    results.reset(new KMeansResult);
    results->centers = centers;
    results->clusterAssignment = clusterAssignment;
    
    return results;
}

/** 
 * Given a collection of cluster centers, this assigns each data point to the closest cluster center
 * and computes the k-means potential for this clustering.
 **/
/*Scalar getKMeansPotential(const vector<Point> &data, const vector<Point> &centers) 
{
    Scalar total = 0;
    vector<int> assignment = assignPointsToClusters(data, centers);
    for (int i = 0; i < (int)data.size(); i++)
        total += distSq(data[i], centers[assignment[i]]);
    return total / data.size();
}*/

void printKMeansResults(KMeansResult* results)
{
    vector<Point>* centers = results->centers.get();
    vector<int>* clusterAssign = results->clusterAssignment.get();
    
    cout << "CENTERS\n";
    for(int i = 0; i < centers->size(); i++)
    {
        cout << (*centers)[i] << "\n";
    }
    cout << "\n\nCLUSTER ASSIGNMENT\n";
    for(int i = 0; i < centers->size(); i++)
    {
        cout << "CLUSTER " << i << "\n";
        for(int j = 0; j < clusterAssign->size(); j++)
        {
            if((*clusterAssign)[j] == i)
            {
                cout << j << " ";
            }
        }
        cout << "\n\n";
    }
}

Scalar euclideanDistance(Point* si, Point* sj)
{
    if(si->getNumDimensions() == sj->getNumDimensions())
    {
        int size = si->getNumDimensions();
        double sum = 0.0;
        for(int i = 0; i < size; i++)
        {
            sum += pow((*si)[i] - (*sj)[i], 2); // Get the sqr
        }
        sum = sqrt(sum);
        return sum;
    }
    else
    {
        return 0.0;
    }
}

void preorder(clustalw::Node* root)
{
    if(root != 0)
    {
        cout << "NODE " << root->seqNum << " HEIGHT " << root->height << "\n"; 
        cout << "LEFT ";
        preorder(root->left);
        cout << "RIGHT ";
        preorder(root->right);
    }
}
