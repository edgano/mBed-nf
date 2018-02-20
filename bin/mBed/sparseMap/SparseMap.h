#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include "pairwise/FullPairwiseAlign.h"
#include "pairwise/FastPairwiseAlign.h"
#include "alignment/Alignment.h"
#include "general/SymMatrix.h"
#include "general/SymMatrix.h"
#include "tree/TreeInterface.h"
#include "tree/UPGMA/RootedGuideTree.h"
#include "tree/UPGMA/RootedClusterTree.h"
#include "sparseMap/Point.h"
#define MAX_COLLISIONS 100
using namespace std;
namespace clustalw
{

class SparseMap
{

    public:
        SparseMap(Alignment* seqs, bool _fast);
	//void doEmbeddingWithSeeds(vector<int> seeds);
        void doEmbedding(string method, int numInputSeeds, int numOutputDims, bool useSeedsOnly,bool findOutliers, vector<int> seeds,string seedSelectionMethod, string seedTreatmentMethod, double distMaskCutoff);

        int getNumSeqs(){return numSeqs;}

        void createRandomGroups();
        void createReferenceGroups(int numInputSeeds, bool findOutliers);
	void createSingletonGroups();
        
	void getSeedSeqsFromUser(vector<int> seeds);
	void getSeedSeqsFromMinimalSampling();
	void getSeedSeqsFromUniformSampling(int numInputSeeds);
	void getSeedSeqsFromRandomSampling(int numInputSeeds);
	void getSeedSeqsFromLengthSampling(int numInputSeeds);
	void getSeedSeqsFromMedoidSampling(int numInputSeeds);


	//void getSeedSeqsFromSingleSampling(int numInputSeeds);
	//void getSeedSeqsFromMedoidsClustering(int numInputSeeds);
        //void selectSeedSeqs(int numInputSeeds, bool findOutliers);


	void reduceNumSeeds(double distMaskCutoff);

	void findOutliersBySeqLengths(double numDeviations);
	void findOutliersBySeqDistances();
	
	void findSeedsUsingPivotObjects();
	void findSeedsUsingPivotGroups();
	
	void findSeedsUsingPivotObjectsAndConverge();
	void findSeedsUsingPivotGroupsAndConverge();


	void findPivotSequence(int seq);
	void findPivotGroup(int seq);
        //void analyseSeedSequences();
        void findSeedSequence(int sequence);
        void identifyOutliers();
         
		
	void findSeedsFromCluster(vector<int> useSeqs);
	void findClustersFromTree(Node* cluster, double cutoff);


	//void findSeedsFromKMedoids(vector<int> useSeqs, int k);


        void findSmallClusters(vector<int> useSeqs);
        void getClusterHeight(Node* cluster);
        void pruneTree(Node* cluster, double cutoff);     
          
        void initialiseCoordinates();
        void generateCoordinates(int numOutputDims,bool useSeedsOnly);
        void reduceDimensionality();
        void greedyResampling(string method);
        void calculateDistanceMatrix();
        void printDistanceMatrixToFile();
        void printCoordinatesToFile();
        void printSequenceOrderToFile();
  
	void addSequenceAsSeed(int sequence, double distMaskCutoff);
	void shellSort(vector<int>* seqList, vector<double>* seqVals);
	void evaluateEmbedding();
        void evaluateEmbedding2();

        auto_ptr<SymMatrix> getDistMat(){return sparseMapDistMat;}
        auto_ptr<vector<Point> > getData(){return data;}

    
    private:
         
	double calculateCorrelation(vector<double> X, vector<double> Y);
	double calculateStress(vector<double> X, vector<double> Y);  
        
        bool vectorContainsElementAlready(vector<int>* group, int element);
	double distMaskCutoff;
	int findLongestSequence();
	int findFurthestObjectFromSequence(int sequence, double distMaskCutoff);
	int findFurthestObjectFromGroup(vector<int> group, double distMaskCutoff);
        int getUniqueNum(vector<int>* group);
		int vectorContainsElementAt(vector<int>* group, int element);
        double treeHeight;
        double getTreeHeight(Node* root);
        double getDistance(int seq1, int seq2);
	double getDistanceFromGroup(int seq, vector<int> group);
        double getMinDistanceInDimension(int seq, int dim);
        int getSeqWithShortestApproxDistance(int seq, int dim);
        double getApproxEuclideanDist(int seq1, int seq2, int dim);
        double getEuclideanDist(Point* vec1, Point* vec2);
        double getEmbeddedDistance(Point* vec1, Point* vec2);
        double getEmbeddedDistance2(Point* vec1, Point* vec2, int k);
        double calcEuclideanDist(vector<double> point1, vector<double> point2);
        double square(double i){return i*i;}
        double dist;
        long double sumDist;
        int numSeqs;
        int beta;
        int kapa;
        int dimensions;
        const float DistNotCalculatedYet;// = -2.0001;
        
        vector<int> useSeqs;
        vector<int> clusterElements;
        vector<int> addGroup1;
        vector<int> addGroup2;
        vector<int> addGroup3;
        vector<int> seqLengths;
        vector<int> seedSeqs;
        vector<string> alignmentNames;
        vector<string> subAlignmentNames;

        vector<vector<int> > dimensionGroups; 
        //vector<string>* sequences;

        auto_ptr<vector<Point> > data;
        auto_ptr<vector<Point> > data_gr;
        SymMatrix distancesFromFast;
        //SymMatrix subDistMat;
        vector<vector<double> > subDistMat;
        auto_ptr<SymMatrix> sparseMapDistMat;
        int counter;
        int counter2; 
        int counter3;
        int counter4; 
        int totalPossible;
        bool fast;
        
	Alignment* alignPtr;

	Alignment seedObj;
        Alignment subAlignObj;

        FastPairwiseAlign fastAlign;
        FullPairwiseAlign fullAlign;

 public:
	int **ppiDistPairCount; /* FS, debug */
	int *piHash; /* FS, debug */
	int iHashKey; /* FS, debug */
	myHash **pprHash;
	int iMaxCollision;
	int ziCollisions[MAX_COLLISIONS];
};

}


