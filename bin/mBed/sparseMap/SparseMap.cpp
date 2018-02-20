#include <cstdlib>
#include <ctime>
#include <math.h>
#include <fstream>
#include <exception>
#include <iostream>
#include <iomanip>
#include <list>
#include <algorithm>
#include "SparseMap.h"
#include "fileInput/FileReader.h"
#include "general/SymMatrix.h"
#include "tree/TreeInterface.h"
#include "tree/UPGMA/RootedGuideTree.h"
#include "tree/UPGMA/RootedClusterTree.h"
namespace clustalw

{
  /* this just reads in the sequences (seqs) and sets up some counters */
  SparseMap::SparseMap(clustalw::Alignment* seqs, bool _fast)
    :dimensions(0),
     DistNotCalculatedYet(-2.0),
  
     counter(0),  // counter  ++ each time the pair-wise alignment function is called;
     counter2(0), // counter2 ++ each time the main getDistance function is called;
     counter3(0), // counter3 ++ each time the distance has already been calculated, and is retrieved from the partial distance matrix (not used anymore)
     counter4(0)  // counter4 ++ each time the two sequences are identical, and dist automatically returned as 0.0.
     
  {
    fast = _fast;
    numSeqs = seqs->getNumSeqs();
    
    alignPtr = seqs;
    kapa = (int)(log2(numSeqs));
    data.reset(new vector<Point>);
    
    int **ppiDistPairCount; /* FS, debug */
    int *piHash; /* FS, debug */
    int iHashKey; /* FS, debug */       
    myHash **pprHash;
    int iMaxCollision;
    int ziCollisions[MAX_COLLISIONS];
  }
  
  void SparseMap::doEmbedding(string method, int numInputSeeds, int numOutputDims, 
			      bool useSeedsOnly, bool findOutliers, vector<int> seeds,
			      string seedSelectionMethod, string seedTreatmentMethod, double distMaskCutoff)
  {
    
  /* Hashing */
  /* a full matrix has N*N entries, 
     because of symmetry on only needs ~ N*(N-1)/2,
     we use (N+1)*N/2 because of the Fortran off-set, 
     as the number of sequences grows the required coverage goes down, 
     the formula 10/sqrt(n) is completely heuristic 
     (guessed sqrt, roughly fitted factor), 
     for small problems there is no point skimping on memory, 
     so use at least BIGISH_PRIME */
  int iFullMatrix = (numSeqs+1)*numSeqs/2;
  double dEstimatedCover = 10/sqrt(numSeqs);
  if (dEstimatedCover > 1.00){
    dEstimatedCover = 1.00;
  }
  int iEstimatedCalls = (int)(dEstimatedCover * iFullMatrix);
#define BIGISH_PRIME 1572869 /* prime half-way between 2^20 and 2^21 */
  iHashKey = (BIGISH_PRIME > iEstimatedCalls) ? BIGISH_PRIME : iEstimatedCalls;

  if (NULL == (pprHash = (myHash **)malloc(iHashKey * sizeof(myHash)))){
    printf("%s:%d: could not malloc %d hash buckets\n", 
	   __FILE__, __LINE__, iHashKey);
    exit(-1);
  }
  memset(pprHash, 0, iHashKey*sizeof(myHash *));

#ifdef DIAGNOSTIC
  iMaxCollision = 0;
  memset(ziCollisions, 0, MAX_COLLISIONS*sizeof(int));
#endif /* diagnostic */


    if(method == "SparseMap")
      {
        createRandomGroups();
      }/*This is the original SparseMap algorithm of Hristescu and Farach-Colton*/
    
    
    if(method == "SeedMap")
      {
	// 1. selection of seed sequences.
	if(seeds.size() > 0) {
	    getSeedSeqsFromUser(seeds);
	   }
	
	else if(seedSeqs.size() == 0)
	  {
	    if(seedSelectionMethod      == "minimal"){getSeedSeqsFromMinimalSampling();}// takes just the longest sequence from the input data
	    else if(seedSelectionMethod == "uniform"){getSeedSeqsFromUniformSampling(numInputSeeds);}// takes a uniformly distributed sampling from the input data
	    else if(seedSelectionMethod == "length") {getSeedSeqsFromLengthSampling(numInputSeeds);}// same as for uniform, but use sequence lengths to order the sequences (very crude)
	    else if(seedSelectionMethod == "random") {getSeedSeqsFromRandomSampling(numInputSeeds);}// takes a randomly distributed sampling from the input data
	    else if(seedSelectionMethod == "medoid") {getSeedSeqsFromMedoidSampling(numInputSeeds);}// takes a randomly distributed sampling from the input data
	    
	    /* 
               TODO:
               ADD FUNCTION FOR PCA/MDS STYLE SEED SELECTION ?
               1. pick a subset of the data (random sampling).
               2. run MDS on this subset, and retrieve eigen vectors.
               3. identify single sequence that most closely mimics the "origin" of each of the top eigen vectors (value closest to 0);
               4. return these sequences as seeds.
	       
               can be run once, or multiple times using different subsets.
	       
	    */
	    
	  }
	
        reduceNumSeeds(distMaskCutoff);

	// 2. treatment of seed sequences. 
	if((seedTreatmentMethod     == "none") || (seedTreatmentMethod =="naive")){cout << "no further seeds will be sought.\n";}
	else if(seedTreatmentMethod == "useLenOutliers"){findOutliersBySeqLengths(3.0);}
	else if(seedTreatmentMethod == "usePivotObjects"){findSeedsUsingPivotObjects();}
	else if(seedTreatmentMethod == "usePivotGroups"){findSeedsUsingPivotGroups();}
	else if(seedTreatmentMethod == "usePivotObjectsAndConverge"){findSeedsUsingPivotObjectsAndConverge();}
	else if(seedTreatmentMethod == "usePivotGroupsAndConverge"){findSeedsUsingPivotGroupsAndConverge();}
	
	createSingletonGroups(); // 3. create the reference groups from the seeds
      }
    
    
    generateCoordinates(numOutputDims, useSeedsOnly);//4. generate the coordinates from the reference groups;
    
    cout << counter  << " times was the pair-wise alignment function called.\n";
    cout << counter2 << " times was the getDistance function called.\n";
    cout << counter4 << " times was the sequence pair identical\n";
    
    cout << "\nFinal Output:[seq/dims]\t" << numSeqs << " x " << (*data)[0].getNumDimensions() << "\n";

#ifdef DIAGNOSTIC
    FILE *pfColl = NULL;
    char zcColl[] = "collisions.dat";
    if (NULL == (pfColl = fopen(zcColl, "w"))){
      printf("%s:%d: could not open file %s for writing\n",
             __FILE__, __LINE__, zcColl);
      exit(-1);
    }
    fprintf(pfColl,
            "#\n"
            "# File: %s\n"
            "#\n"
            "# Note: generated by %s:%d\n"
            "#       %d hash buckets\n"
            "#\n"
            , zcColl, __FILE__, __LINE__, iHashKey);
    for (int i = 0; i < MAX_COLLISIONS; i++){
      fprintf(pfColl, "%d\t%d\n", i, ziCollisions[i]);
    }
#endif /* diagnostic */


    /* free up (hash) resources */
    if (NULL != pprHash){
      for (int i = 0; i < iHashKey; i++){
	if (NULL != pprHash[i]){
	  free(pprHash[i]);
	  pprHash[i] = NULL;
	}
      }
      free(pprHash);
      pprHash = NULL;
    }

	

  }
  
  void SparseMap::createRandomGroups()
  {    
    /* 
       The original SparseMap algorithm uses a collection of log2**2 random groups
       to act as reference groups (stored in the vector "dimensionGroups") 
    */
    try
      {
        dimensionGroups.clear();
        srand((unsigned)time(0));
        int num;
        bool exists;
        int groupSize = 2;
        bool foundUnique;
	
        cout << "numSeqs = " << numSeqs << "=> kapa = " << kapa << " => numDims = " << kapa * kapa << "\n";    
        for(int i = 0; i < kapa; i++) // kapa rows
	  {
            cout << "Generating row " << i + 1 << " of " << kapa << "\n";
            for(int j = 0; j < kapa; j++) // kapa elements in each row
	      {
                vector<int> group(groupSize, -1);
                for(int l = 0; l < group.size(); l++)
		  {
                    num = getUniqueNum(&group);   
                    group[l] = num;
		  }
                dimensionGroups.push_back(group);            
	      }
            groupSize *= 2;
            if(groupSize == numSeqs)
	      {
                break;
	      }
	  }
      }
    catch(const exception &ex)
      {
        cerr << "Exception in createRandomGroups function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  } 
  int SparseMap::getUniqueNum(vector<int>* group)
  {
    try
      {
	int num;
	bool vectorContains = true;
	while(vectorContains)
	  {
	    num = (rand() % numSeqs);
	    vectorContains = vectorContainsElementAlready(group, num);
	  }
	return num;
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getUniqueNum function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }    
  }
  void SparseMap::createSingletonGroups()
  {    
    try
      {
	
        dimensionGroups.clear(); 
	useSeqs.clear(); 
	addGroup1.clear();         
        
	for(int i=0; i < seedSeqs.size(); i++)
	  {
	    addGroup1.push_back(seedSeqs[i]);
	    dimensionGroups.push_back(addGroup1);
	    addGroup1.clear();
	  }     
	
	seedSeqs.clear();
	
	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in createSingletonGroups function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }        
  } 
  
  
  /******
	 
	 the following functions are the various seed selection methods.
	 Currently the options are to:
	 
	 1. pick just the longest sequence (really recommend you follow it up with some form of seedTreatment.
	 2. pick pow(log2N,2) sequences uniformly from the input set
	 3. pick pow(log2N,2) sequences randomly from the input set
	 4. pick pow(log2N,2) sequences uniformly from the input set (sorted by sequence length)
	 5. let the user decide which sequences to use (NOT WORKING SINCE ANSCHLUSS)
	 
	 NOTE: since ANSCHLUSS, sequences are no longer automatically 
	 sorted by sequence length. This means that, for example, the 
	 order of the input sequences not has an effect on the sequences
         that are selected by these sampling methods. It's not a problem,
	 as I have an option in place to fix these sequences in place (the 
	 getSeedSeqsFromLengthSampling function), if someone wants to use it. 
  ******/
  
  
  void SparseMap::getSeedSeqsFromUser(vector<int> seeds)
  {

    /* 
       this takes the input from a seed file (which has been parsed in main.cpp).
       NOTE: Since ANSCHLUSS, this is no longer possible (for now). 
    */
    try
      {
	seedSeqs.clear();
	//---
	for(int i = 0; i < seeds.size(); i++)
	  {
	    seedSeqs.push_back(seeds[i]);
	  }
	//----
	seeds.clear();
	cout<< "SparseMap numSeeds: " << seedSeqs.size() << "\n\n";
      }
    
    catch(const exception &ex)
      {
	cerr << "Exception in getSeedSeqsFromUser function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }
  }
  void SparseMap::getSeedSeqsFromMinimalSampling()
  {
    // this selection process is simple (Longest sequence only)
    try
      {
        seedSeqs.clear();
        int longestSeq = findLongestSequence();
	seedSeqs.push_back(longestSeq);
	
        cout<< "using " << seedSeqs.size() << " sequence from minimal sampling (longest sequence only).\n";
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getSeedSeqsFromMinimalSampling function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  } 
  void SparseMap::getSeedSeqsFromUniformSampling(int numInputSeeds)
  {
    // this selection process is simple (uniform sampling idea)
    try
      {
        seedSeqs.clear();
	// if the number of input seeds is not specified, 
	// get pow(kapa,2) seeds from data (kapa = log2N)
	
	int numSeeds = 0;
	if(numInputSeeds == -1){numSeeds = kapa * kapa;}
	else if(numInputSeeds >=0){numSeeds = numInputSeeds;}
	
	double numBisections = log2(numSeeds);
	
	while(pow(numBisections,2) < numSeeds){numBisections++;}
	
	int numScans = (int)(numBisections + 1);
	
	//cout << numSeeds << " seeds requires " << numScans << " scans of the input data\n";
	
	for(int scan = 1; scan < numScans+1; scan++)
	  {
	    int p = (int)(pow(2.0,scan));
	    
	    for(int j = 1; j < p; j+=2)
	      {
		int seed = j * numSeqs / p;
		
		if(seedSeqs.size() < numSeeds && vectorContainsElementAlready(&seedSeqs,seed) == false)
		  {
		    seedSeqs.push_back(seed);
		  }
	      }
	  }
        cout<< "using " << seedSeqs.size() << " sequences from uniform sampling\n";
	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getSeedSeqsFromUniformSampling function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  }  
  void SparseMap::getSeedSeqsFromRandomSampling(int numInputSeeds)
  {
    // this selection process is simple (random sampling idea)
    try
      {
        seedSeqs.clear();
        srand((unsigned)time(0));
        int num;
        bool exists;
        bool foundUnique;
        
	// if the number of input seeds is not specified, 
	// get pow(kapa,2) seeds from input (kapa = log2N)
	
	int numSeeds = 0;
	if(numInputSeeds == -1)   {numSeeds = kapa * kapa;}
	else if(numInputSeeds >=0){numSeeds = numInputSeeds;}
	if(numSeeds > numSeqs)    {numSeeds = numSeqs;}
	
	cout << "looking for " << numSeeds << " seed sequences\n";
	
	for(int i = 0; i < numSeeds; i++)
	  {
	    int seed = getUniqueNum(&seedSeqs);
	    seedSeqs.push_back(seed);
	  }
        cout<< "\nusing " << seedSeqs.size() << " sequences from random sampling\n";
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getSeedSeqsFromRandomSampling function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  } 
  void SparseMap::getSeedSeqsFromLengthSampling(int numInputSeeds)
  {
    // this selection process is simple (uniform sampling idea, but sorted by sequence lengths)
    try
      {
        seedSeqs.clear();
	
	// need a list of the input sequences, sorted by sequence length, so use the shellSort function. 
	// shellSort works with doubles, so make each length a float. 
	
        double Len = 0;   
	
        vector<int> seqList;
        vector<double> seqLens;
        
	for(int Seq=0; Seq < numSeqs; Seq++)
	  {
	    Len = (double)(alignPtr->getSeqLength(Seq + 1));
	    seqList.push_back(Seq);
	    seqLens.push_back(Len);	
	  }
	
	shellSort(&seqList, &seqLens);
	
	// if the number of input seeds is not specified, 
	// get pow(kapa,2) seeds from data (kapa = log2N)
	
	int numSeeds = 0;
	if(numInputSeeds == -1)   {numSeeds = kapa * kapa;}
	else if(numInputSeeds >=0){numSeeds = numInputSeeds;}
	
	double numBisections = log2(numSeeds);
	while(pow(numBisections,2) < numSeeds){numBisections++;}
	
	
	int numScans = (int)(numBisections + 1);
	
	//cout << numSeeds << " seeds requires " << numScans << " scans of the input data\n";
	
	
       	for(int scan = 1; scan < numScans+1; scan++)
	  {
	    int p = (int)(pow(2.0,scan));
	    
	    for(int j = 1; j < p; j+=2)
	      {
		int seed = j * numSeqs / p;
		
		if(seedSeqs.size() < numSeeds && vectorContainsElementAlready(&seedSeqs,seed) == false)
		  {
		    seedSeqs.push_back(seqList[seed]);
		  }
	      }
	  }
        seqLens.clear();     
        seqList.clear();
        cout<< "using " << seedSeqs.size() << " sequences from length sampling\n";
	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getSeedSeqsFromLengthSampling function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  }  
  void SparseMap::getSeedSeqsFromMedoidSampling(int numInputSeeds)
  {
    try
      {
	
	/*	
                PLEASE NOTE: THIS CODE IS TERRIBLE!! I'll need to work on it a lot...
		
		
		This function takes a sample from the input data and quickly finds 
		potential seeds from them using a 'k-means like k-medoids algorithm' 
		(Hae-Sang Park, Jong-Seok Lee and Chi-Hyuck Jun, 2006).
		
                The centroids from this method (or at least the first part) are used as seeds.
		
		If this is run several times using different samplings, should get a fairly good coverage. 
	*/ 
	
	seedSeqs.clear();  // seedSeqs holds the final seed selection. 
	addGroup1.clear(); // addGroup1 holds the sequences to be selected at each iteration.
	addGroup2.clear(); // addGroup2 holds the sequences that have been selected overall. 
	
	vector<int> tempGroup; // tempGroup contains the same elements as addGroup1 - but it can be sorted, resorted, etc...
	vector<int> centroids;
	
	
	double dist, centroidDistance, globalSumOfDistances;
	int seq1, seq2, seed, centroid, ind1, ind2;

	
	for(int iter = 0; iter < kapa; iter++)
	  {
	    addGroup1.clear(); centroids.clear(); tempGroup.clear();
	    
	    
	    centroidDistance = 2.0;
	    globalSumOfDistances = 0.0;
	    
	    
	    /*	
		STEP1:		(selection of initial medoids)
		
		step_1_1:
		
		generate a full d(i,j) distance matrix. 
		For this we need some sequences. May as well select some random seqs each time. 
		
	    */
	    
	    int nSeqs = 0;
	    if(numInputSeeds == -1)    {nSeqs = kapa * kapa;}
	    if(numInputSeeds >=0)      {nSeqs = numInputSeeds;}
	    if(nSeqs > numSeqs)        {nSeqs = numSeqs;}

	    cout << "looking for " << nSeqs << " sequences\n";
	    
	    for(int i = 0; i < nSeqs; i++)
	      {
		seed = getUniqueNum(&addGroup1);
		if(vectorContainsElementAlready(&addGroup2, seed) == false)
		  {
		    addGroup1.push_back(seed); 
		    addGroup2.push_back(seed);
		    tempGroup.push_back(seed);
		  }
	      }
	    
	    
	    nSeqs = addGroup1.size(); cout<< "\nSelected " << nSeqs << " sequences as input to medoid sampling\n";
	    int numCentroids = (int)(pow(nSeqs,0.5));
	    
	    
	    
	    double dij = 0.0; double dijSum = 0.0;
	    double dil = 0.0; double dilSum = 0.0;
	    double pij = 0.0; double pijSum = 0.0;
	    
	    
	    vector<vector<double> > dijDistMat; for(int i = 0; i < nSeqs; i++){dijDistMat.push_back(vector<double>(nSeqs,2.0));}

	    for(int i = 0; i < nSeqs; i++)
	      {
		for(int j = 0; j <= i; j++)
		  {
		    dij = getDistance(addGroup1[i], addGroup1[j]); dijDistMat[i][j] = dij;  dijDistMat[j][i] = dij;
		  }  
	      }
	    
	    
	    /*	step_1_2
		
		calculate p(i,j) =[d(i,j) / sum (d(i,l))] to make an initial guess at cluster centres. I take l
		to be any third object in the input data (i.e. we could use j, but wanna keep the variables separate) 
	    */
	    
	    vector<vector<double> > pijDistMat; for(int i = 0; i < nSeqs; i++){pijDistMat.push_back(vector<double>(nSeqs,2.0));}
	    
	    for(int i = 0; i < nSeqs; i++)
	      {
		dilSum = 0.0;
		for(int l = 0; l < nSeqs; l++)
		  {
		    dil = dijDistMat[i][l];
		    dilSum += dil;
		  }
		for(int j = 0; j < nSeqs; j++)
		  {
		    dij = dijDistMat[i][j];
		    pij =  dij / dilSum;
		    pijDistMat[i][j] = pij;
		    pijSum += pij;
		  }  
	      }
	    
	    /* step_1_3:	
	       calculate sum(pij)) for all j objects, and sort in ascending order. 
	       
	       Select first k objects having minimal values of p(i,j) to act as initial medoids.
	    */
	    
	    vector<double> pijVals;
	    
	    for(int i = 0; i < nSeqs; i++)
	      {
		pijSum = 0.0;
		for(int j = 0; j < nSeqs; j++)
		  {
		    pijSum += pijDistMat[i][j];
		  }
		pijVals.push_back(pijSum);
	      }
	    
	    
	    
	    shellSort(&tempGroup, &pijVals);// I want to keep addGroup1 as it is - no sorting please - use tempGroup instead!! 
	    
	    
	    
	    
	    /*	step_1_4: 
		assign each object to nearest centroid.
		For every object i, calculate distance to every centroid j. 
		can do this in a distance matrix. For each row, find minimal value. 
		
		The initial choice of centroids is the first k elements in tempGroup 
		By default, k is set at sqrt(nSeqs)... this is equivalent to kapa (if kapa is used above) 
		
		
	    */	
		
	    // dicDistMat is the distance matrix for distances between every object i and every centroid c
	    vector<vector<double> > dicDistMat; for(int i = 0; i < nSeqs; i++){dicDistMat.push_back(vector<double>(numCentroids,2.0));}
	    
	    
	    vector<int> centroidAssignments;
	    
	    for(int i = 0; i < nSeqs; i++)
	      {
		seq1 = addGroup1[i]; 
		ind1 = vectorContainsElementAt(&addGroup1, seq1);
		
		centroidDistance = 2.0;
		
		for(int j = 0; j < numCentroids; j++)
		  {
		    seq2 = tempGroup[j];
		    ind2 = vectorContainsElementAt(&addGroup1, seq2);
		    
		    dist = dijDistMat[ind1][ind2]; 
		    dicDistMat[i][j] = dist;
		    
		    // cout << "looking at indices " << seq1 << "," << seq2 << " ( distance = " << dist << ")\n";
		    if(dist < centroidDistance)
		      {
			centroid = seq2; 
			centroidDistance = dist;
		      }
		  }
		centroidAssignments.push_back(centroid);
		if(vectorContainsElementAlready(&centroids, centroid) == false){centroids.push_back(centroid);}
		
		cout << " seq " << seq1 << ": assigned to " << centroid << " (distance = " << centroidDistance << ")\n";
	      }
	    
	    
	    /*	step_1_5	
		
		Calculate current optimal value (sum of distances from all objects to the nearest centroid)
	    */
	    
	    globalSumOfDistances = 0.0;
	    for(int i = 0; i < centroids.size(); i++)
	      {
		
		for(int j = 0; j < nSeqs; j++)
		  {
		    seq1 = centroids[i];
		    seq2 = addGroup1[j];
		    
		    dist = dijDistMat[vectorContainsElementAt(&addGroup1, seq1)][vectorContainsElementAt(&addGroup1, seq2)];
		    
		    globalSumOfDistances += dist;
		    
		  }
		
	      }
	    
	    cout << "globalSumOfDistances = " << globalSumOfDistances << "\n";
	    
	    
	    
	    /*	STEP2:		(Find new medoids)
		
		Replace current medoid in each cluster with object
		that minimises distance to all other objects in cluster
	    */
	    
	    for(int i = 0; i < nSeqs; i++){dicDistMat[i] = vector<double>(numCentroids,2.0);}
	    
	    centroidDistance = 2.0;
	    dist = 0.0;
	    centroid = 0;
	    
	    
	    
	    for(int j = 0; j < centroids.size(); j++)
	      {  
		vector<int> clusterMembership;
		
		centroid = centroids[j];
		cout << centroid << ":\t(";
		for(int i = 0; i < nSeqs; i++)
		  {
		    if(centroid == centroidAssignments[i])
		      {
			clusterMembership.push_back(addGroup1[i]); 
			cout << addGroup1[i] << ", ";
		      }
		  }
		cout <<")\n";
		double minRowSum = 100.0;
		int newCentroid = centroids[j];
		
		int cSeqs = clusterMembership.size();
		vector<vector<double> > clusDistMat;for(int i = 0; i < cSeqs; i++){clusDistMat.push_back(vector<double>(cSeqs,2.0));}
		
		double currentSumOfDistances = 0.0;
		for(int i = 0; i < cSeqs; i++)
		  {
		    seq1 = vectorContainsElementAt(&addGroup1, centroid);
		    seq2 = vectorContainsElementAt(&addGroup1, clusterMembership[i]);
		    
		    currentSumOfDistances += dijDistMat[seq1][seq2];
		  }
		cout << centroid << " is the current centroid (rowSum = " << currentSumOfDistances <<")\n";
		
		for(int i = 0; i < cSeqs; i++)
		  {
		    
		    seq1 = clusterMembership[i];
		    
		    double rowSum = 0.0;
		    for(int j = 0; j < cSeqs; j++)
		      {
			seq2 = clusterMembership[j];
			dist = dijDistMat[vectorContainsElementAt(&addGroup1, seq1)][vectorContainsElementAt(&addGroup1, seq2)];
			clusDistMat[i][j] = dist;
			rowSum += dist;
		      }
		    cout << seq1 << " " << seq2 << " " << rowSum <<"\n";
		    
		    if(rowSum <= currentSumOfDistances)
		      {
			currentSumOfDistances = rowSum; newCentroid = seq1;
		      }
		  }
		//cout << newCentroid << " is selected as new centroid (rowSum = " << currentSumOfDistances <<" \n";
		centroids[j] = newCentroid; 
	      }
	    cout <<"ading the following sequences as seeds: ";
	    for(int a = 0; a < centroids.size(); a++)
	      { 
		cout << centroids[a] <<" ";seedSeqs.push_back(centroids[a]);
	      }
	    cout <<"\n";
	    
	    ////	STEP3:		(New Assignment)
	    
	    ////	step_3_1	Assign each object to nearest medoid
	    
	    ////	step_3_2	Calculate new optimal value, the sum of distance from all objects to their new medoids.
	    ////		      	If new optimal value is equal to previous optimal value, stop algorithm.
	    //                          otherwise, go back to step 2
	    dijDistMat.clear();
	    pijDistMat.clear();
	    dicDistMat.clear();
	    addGroup1.clear(); pijVals.clear();
	    centroidAssignments.clear();
	    
	  }
      }
    catch(const exception &ex)
      {
        cerr << "Exception in getSeedSeqsFromMedoidSampling function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }
  
  /******
	 the following functions are the various seed treatment methods.
	 Currently the options are to:
	 
	 1. reduce the number of seeds in the list by comparing them all and rejecting the similar ones (this is done by default)
	 2. finding outliers by sequence lengths
	 3. finding seeds using pivot Objects, or finding seeds using pivot Groups. These are explained further below. 
	 4. finding seeds by (3), but allowing them to run to completion rather than setting a predefined cutoff. 
  ******/
  
  void SparseMap::reduceNumSeeds(double distMaskCutoff)
  {
       /*
      This function compares all the seeds that have been chosen. If two seeds are very similar
       (below the distMaskCutoff value), the longer one is retained, the other one is discarded.

      Right now the cutoff value has to be specified on the command line, but I want an automated
      way of doing this too. Something along the lines of:
      1. get distances between all seeds.
      2. use these distances to specify a particular cutoff value. Perhaps by finding the
      mean, std dev, and using those, or percentiles, etc.
      */

      try
      {
        cout << "Restricting seed selection using a distance cutoff of " << distMaskCutoff << ".\n";
	cout << "if two sequences are closer than this cutoff, the shorter one is discarded \n";

	int nSeeds = seedSeqs.size();
	addGroup1.clear(); 
	for(int i  = 0; i < nSeeds; i++){addGroup1.push_back(seedSeqs[i]);}
	seedSeqs.clear();
	
	double dij;
	int seq1, seq2;
	
        vector<vector<double> > dijDistMat;
        for(int i = 0; i < nSeeds; i++){dijDistMat.push_back(vector<double>(nSeeds,2.0));}
	
        for(int i = 0; i < nSeeds; i++)
	  {
	    
	    for(int j = 0; j <= i; j++)
	      {
		seq1 = addGroup1[i];
		seq2 = addGroup1[j];
		
		if(seq1 == seq2){dij = 0.0;}
		
		else if(seq1!=seq2)
		  {
		    dij = getDistance(seq1, seq2); 
		    if(dij > distMaskCutoff){dij = dij;}
		    else if(dij <=distMaskCutoff){dij = 2.0;}
		    
		  }
		
		dijDistMat[i][j] = dij;  dijDistMat[j][i] = dij;
	      }
	    
	  }
	for(int i = 0; i < nSeeds; i++)
	  {
	    for(int j = 0; j <= i; j++)
	      {
		seq1 = addGroup1[i];
		seq2 = addGroup1[j];
		
		dij = dijDistMat[i][j];
		
		if(dij == 2.0)
		  {
		    int len1 = alignPtr->getSeqLength(seq1 + 1); 
		    int len2 = alignPtr->getSeqLength(seq2 + 1);
		    if(len1 >= len2){useSeqs.push_back(seq2);}
		    if(len1 <  len2){useSeqs.push_back(seq1);}
		  }
		
	      }  
	  }
	
	for(int i = 0; i < nSeeds; i++)
	  {
	    seq1 = addGroup1[i]; 
	    if(vectorContainsElementAlready(&useSeqs, seq1) == false){seedSeqs.push_back(seq1);}
	  }		
	addGroup1.clear();dijDistMat.clear(); 
	
      }
      catch(const exception &ex)
	{
	  cerr << "Exception in reduceNumSeeds function\n";
	  cerr << ex.what() << endl;
	  cerr << "Terminating program. Cannot continue\n";
	  exit(1);
	}        
  }   
  void SparseMap::findOutliersBySeqLengths(double numDeviations)
  {
    /*
      This function uses the lengths of the sequences to find potential outliers
      using simple things like average length and standard deviations. This is
      probably very weak from a mathematical/logical point of view, though.
    */
    try
      {
        cout << "selecting outliers from length distribution\n";
	
	// first off, calculate the average and standard deviation of the sequence lengths.
	unsigned long int sumLen = 0; // for large N, this is gonna get big, need the extra memory of an unsigned long int
	
	int seqLen = 0;
	for(int i = 0; i < numSeqs; i++)
	  {
	    sumLen += alignPtr->getSeqLength(i + 1);       
	  }
	float aveLen = sumLen / numSeqs;
	
	sumLen = 0;
	
	for(int i = 0; i < numSeqs; i++)
	  {
	    sumLen += (int)(pow((alignPtr->getSeqLength(i + 1) - aveLen),2));
	  }
	float StdDev = sqrt(sumLen / (numSeqs - 1));
	
	
	/* 
	   The boundaries listed below are near the limits of normal distribution.
	   Therefore, if we were to assume that the sequences are normally distributed (e.g.): 
	   
	   if numDeviations = 2.5 : ~ 98% of the sequences should be within this boundary
	   if numDeviations = 3.0 : ~ 99.7% of the sequences should be within this boundary.
	   
	   
	   If more than one sequence has the same length, we just take the first one
	   that comes along. This was added because of the huge number of sequences
	   that I was finding on big test cases. This analysis is based on length only, 
	   and as such, sequences of equal length should be considered identical. 
	   
	   Another idea (not implemented yet) is that for each length outside the 
	   boundary, grab all sequences, and find a more suitable representative 
	   sequence. Gonna test this by clustering all outlier sequences, and adding just 
	   a single representative from each of the groups found (longest seq).
	*/
	
	
	double lowerBoundary = aveLen - (numDeviations * StdDev);
	double upperBoundary = aveLen + (numDeviations * StdDev);
	
	int Outlier = 0;
	
	vector<int> usedLens;
	
	for(int i = 0; i < numSeqs; i++)
	  {
	    seqLen = alignPtr->getSeqLength(i + 1);
	    
	    // don't want EVERY sequence outside the boundaries; 
	    // just one for each length will do.
	    
	    if(vectorContainsElementAlready(&usedLens, seqLen) ==false) 
	      {
		if(	(seqLen < lowerBoundary) || 
			(seqLen > upperBoundary) &&
			(vectorContainsElementAlready(&seedSeqs, i) == false))
		  
		  {
		    Outlier++;
		    addSequenceAsSeed(i, distMaskCutoff);
		    usedLens.push_back(seqLen);
		  }
	      }
	  }
	usedLens.clear();      
      }       
    catch(const exception &ex)
      {
	cerr << "Exception in findOutliersBySequenceLength function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }        
  }                 
  void SparseMap::findSeedsUsingPivotObjects()
  {
    /*
      pivotObjects are single sequences in the data set.
      This works in a similar way to the choose-distant-objects function in FastMap. Take the
      first sequence, and find the sequence that's furthest away; then find the sequence that's
      furthest away from that sequence
    */
    try
      {
	
	
	int numSearches = seedSeqs.size();
	
	for(int i = 0; i < numSearches; i++)
	  {
	    int seq1 = seedSeqs[i];
	    if(vectorContainsElementAlready(&useSeqs, seq1) == false)
	      {
		cout << "\nAnalysing seed sequence " << alignPtr->getName(seq1 + 1) << " (" << i+1 << " / " << numSearches << ")...";
		findPivotSequence(seq1);
	      }		
	  }	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in findSeedsUsingPivotObjects function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      } 
  }
  void SparseMap::findSeedsUsingPivotGroups()
  {
    /*
      pivotGroups are groups of sequences in the data set.
      This works in a similar way to the centroid identification stage of k-medoids. Take the
      first sequence, and find the sequence that's furthest away; then find the sequence that's
      furthest away from BOTH the first two sequences, etc....
    */
    
    try
      {
	
	
	int numSearches = seedSeqs.size();	
	for(int i = 0; i < numSearches; i++)
	  {
	    int seq1 = seedSeqs[i];
	    if(vectorContainsElementAlready(&useSeqs, seq1) == false){
	      cout << "Analysing seed sequence " << alignPtr->getName(seq1 + 1) << " (" << i+1 << " / " << numSearches << ")...";
	      findPivotGroup(seq1);}
	    
	  }	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in findSeedsUsingPivotGroups function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      } 
  }
  void SparseMap::findSeedsUsingPivotObjectsAndConverge()
  {
    /*
      
      This function works in the same way as usePivotObjects, above.
      The only difference is that each of the sequences identified is subjected
      to the same analysis as the original seeds.
      
    */
    
    try
      {
	for(int i = 0; i < seedSeqs.size(); i++)
	  {
	    int seq1 = seedSeqs[i];
	    if(vectorContainsElementAlready(&useSeqs, seq1) == false){
	      cout << "Analysing seed sequence " << alignPtr->getName(seq1 + 1) << " (" << i+1 << " / " << seedSeqs.size() << ")...";
	      findPivotSequence(seq1);}		
	  }	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in findSeedsUsingPivotObjectsAndConverge function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      } 
  }
  void SparseMap::findSeedsUsingPivotGroupsAndConverge()
  {
    /*
      
      This function works in the same way as usePivotGroups, above.
      The only difference is that each of the sequences identified is subjected
      to the same analysis as the original seeds.
      
    */
    
    try
      {
	for(int i = 0; i < seedSeqs.size(); i++)
	  {
	    int seq1 = seedSeqs[i];
	    if(vectorContainsElementAlready(&useSeqs, seq1) == false)
	      {
		cout << "\n analysing seed sequence " << alignPtr->getName(seq1 + 1) << " (" << i+1 << " / " << seedSeqs.size() << ")...";
		findPivotGroup(seq1);
	      }
	  }	
      }
    
    catch(const exception &ex)
      {
	cerr << "Exception in findSeedsUsingPivotGroupsAndConverge function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      } 
  }
  void SparseMap::findPivotSequence(int sequence)
  {
    /*
      this function analyses a sequence. It does this by:
      1. calculating distances between the sequence and all other sequences.
      The sequence that is furthest away is called a pivot object (pivotA).
      2. step 1 is repeated with pivotA, to find pivotB.
      //3. step 1 is repeated with pivotB, to find pivotC.
      
      
      |      *       Take note that the two most distant points from a given sequence are not
      |    /   \     necessarily close together, even if that sequence is a pivot.
      |   /     \    Think of an equilateral triangle: from the top corner,
      |  /       \   the distances to other two corners are identical.
      | *---------*

    */


    try
      {
	
	
        /*
	  1. if the seed sequence has not previously been identified as a pivot,
	  we analyse the seed sequence, and find the most distant point, pivotA
	  The seed should always be analysed; the later "if" statements should ensure it!
        */
	
	//if(vectorContainsElementAlready(&useSeqs, sequence) == true) { cout << " has already been analysed. Stop.\n";}
	
	
	if(vectorContainsElementAlready(&useSeqs, sequence) == false)
	  {
	    //cout << "...analysing\n";
	    int pivotA = findFurthestObjectFromSequence(sequence, distMaskCutoff);
	    //cout<< "pivotA = " << alignPtr->getName(pivotA + 1); 
            
	    useSeqs.push_back(sequence);
            
	    
	    
            /* 2. if pivotA has previously been analysed, or is queued, stop.*/
	    
	    /*
	      if(vectorContainsElementAlready(&seedSeqs, pivotA) == true)
	      {
	      cout << " is also a seed sequence";
	      if(vectorContainsElementAlready(&useSeqs, pivotA) == true)
	      {
	      cout << " and has already been analysed. Stop.\n";
	      }
	      else
	      {
	      cout << " and is queued to be analysed. Stop.\n";
	      }
	      }
	      
	      
	      else if((vectorContainsElementAlready(&useSeqs, pivotA) == true) &&
	      (vectorContainsElementAlready(&seedSeqs,pivotA) == false))
	      {
	      cout << " has already been analysed as a pivot. Stop.\n";
	      }
	    */
	    
	    if ((vectorContainsElementAlready(&useSeqs, pivotA) == false) &&
		(vectorContainsElementAlready(&seedSeqs,pivotA) == false))
	      {
		//cout << "...analysing\n";
		int pivotB = findFurthestObjectFromSequence(pivotA, distMaskCutoff); 
		useSeqs.push_back(pivotA);
		
		// 3. if pivotB has previously been analysed, or is queued, stop. 
                
		//cout << "pivotB = " <<alignPtr->getName(pivotB + 1);
		
		/*
		  if(vectorContainsElementAlready(&seedSeqs, pivotB) == true)
		  {
		  cout << " is also a seed sequence";
		    if(vectorContainsElementAlready(&useSeqs, pivotB) == true)
		    {
		    cout << " and has been analysed already. Stop.\n";
		    }
		    else
		    {
		    cout << " and is queued to be analysed. Stop.\n";
		    }
		    }
		    else if((vectorContainsElementAlready(&useSeqs, pivotB) == true) &&
		    (vectorContainsElementAlready(&seedSeqs, pivotB) == false))
		    {
		    cout << " has already been analysed as a pivot. Stop.\n";
		    }
		*/
		if((vectorContainsElementAlready(&useSeqs, pivotB) == true) &&
		   (vectorContainsElementAlready(&seedSeqs, pivotB) == false))
		  {
		    //seedSeqs.push_back(pivotB);
		    addSequenceAsSeed(pivotB, distMaskCutoff);
		  }
		
	      }
	  }     
      }
    catch(const exception &ex)
      {
	cerr << "Exception in findPivotSequence function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }        
  }     
  int SparseMap::findFurthestObjectFromSequence(int Seq1, double distMaskCutoff)
  {
   /*
    this algorithm is inspired by FastMap. FastMap is a heuristic of PCA, which computes
    eigenvectors from an input set. FastMap tries to mimic this by finding the two
    objects in the dataset with the greatest distance between them. First it chooses
    a random object, and finds the point that is furthest away from it. this is called
    pivot A. Then, from pivot A, it finds the furthest object, which is then called
    pivot B. This choice of pivots depends on the choice of initial random object,
    and so it's recommended to be run several times from different start points.
  */
    try
      {                     
        double dist = 0.0;   
        vector<int> seqList; // this will hold the list of sequences
        vector<double> seqDist; // this will hold a list of the inter-sequence distances
        
        for(int Seq2=0; Seq2 < numSeqs; Seq2++)
	  {   
	    if(vectorContainsElementAlready(&useSeqs, Seq2) == false)
	      {
		dist = getDistance(Seq1, Seq2);
		seqList.push_back(Seq2);
		seqDist.push_back(dist);
		
		/* using the distance masking idea to reduce the number of potential seeds */
		if(dist < distMaskCutoff)
		  {
		    int len1 = alignPtr->getSeqLength(Seq1 + 1); 
		    int len2 = alignPtr->getSeqLength(Seq2 + 1);
		    if(len1 >= len2){useSeqs.push_back(Seq2);}
		    if(len1 <  len2){useSeqs.push_back(Seq1);}
		  }
	      }
	  }  
	shellSort(&seqList, &seqDist);
	int furthestObject  = seqList[seqList.size() - 1]; 
        seqDist.clear();     
        seqList.clear();
        return furthestObject;    
      }
    
    catch(const exception &ex)
      {
        cerr << "Exception in findFurthestObjectFromSequence function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  }   
  void SparseMap::findPivotGroup(int sequence)
  {
    
    try
      {
	/* 
	   1. find the sequence that is at greatest distance from the initial sequence, SEED. This is called pivotA. 
	   2. find the sequence that is at greatest distance from the group (SEED,pivotA). This is called pivotB. 
	   3. find the sequence that is at greatest distance from the group (SEED,pivotA,pivotB)... etc etc etc
	   
	   when one of the identified sequences has already been identified as a pivot sequence or seed, stop. 
	   when the group size reaches a maximal limit (for example, log2(N)), stop.
	*/
	
	// 1. if the seed sequence has not previously been identified as a pivot,
	// we analyse the seed sequence, and find the most distant point, pivotA
	// The seed should always be analysed; the later "if" statements should ensure it!
        
	
	
	/*
	  if(vectorContainsElementAlready(&useSeqs, sequence) == true)
	  {
	  cout << " has already been analysed. Stop.\n";
	  }
	*/
	
	if(vectorContainsElementAlready(&useSeqs, sequence) == false)
	  {
	    //cout << "...analysing\n";
	    vector<int> pivotGroup;
	    int pivotSeq;
	    int maxGroupSize = kapa;
	    int oldSeedNums = seedSeqs.size();
	    
	    
	    pivotGroup.push_back(sequence); 
	    
	    for(int i = 0; i < maxGroupSize; i++)
	      {
		pivotSeq = findFurthestObjectFromGroup(pivotGroup, distMaskCutoff);
		
		
		/*
		  if(vectorContainsElementAlready(&seedSeqs, pivotSeq) == true)
		  { 
		  cout << " is also a seed sequence";
		  
		  if(vectorContainsElementAlready(&useSeqs,  pivotSeq) == true)  
		  { 
		  cout << " and has already been analysed. Stop.\n";
		  }
		  else
		  { 
		  cout << " and is queued to be analysed. Stop.\n";
		  }
		  }
		  
		  

		  if((vectorContainsElementAlready(&useSeqs, pivotSeq) == true) &&				
		  (vectorContainsElementAlready(&seedSeqs,pivotSeq) == false))
		  
		  {
		  cout << " has already been analysed. Stop.\n";
		  }
		*/
		
		if((vectorContainsElementAlready(&useSeqs, pivotSeq)   == false) &&
		   (vectorContainsElementAlready(&seedSeqs,pivotSeq)   == false) &&
		   (vectorContainsElementAlready(&pivotGroup,pivotSeq) == false) )
		  {
		    pivotGroup.push_back(pivotSeq);
		    addSequenceAsSeed(pivotSeq, distMaskCutoff);
		
		  }
	      }
	    useSeqs.push_back(sequence);
	    cout << "added " << seedSeqs.size() - oldSeedNums << " new seed sequences to analysis\n";
	    
	    pivotGroup.clear();
	  }
	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in findPivotGroup function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }        
  }     
  int SparseMap::findFurthestObjectFromGroup(vector<int> group, double distMaskCutoff)
  {    
    try
      {                     
	/* 
           I'm not sure if this is the correct way to do this, but it is A way for now;
	   To find the sequence with the greatest distance from all other sequences in
	   a group, I'm copying the Euclidean distance function. This is wrong cos the 
	   distances are not Euclidean, but other than simply using the sum of all the
	   distances, I don't have many other ideas...
	*/
	
        double dist = 0.0;
	double squareDist = 0.0;
	
        vector<int>    seqList;
        vector<double> seqDist;
        
        for(int seq1=0; seq1 < numSeqs; seq1++)
	  {
	    squareDist = 0.0;
	    for(int i=0; i < group.size(); i++)
	      {
		int seq2 = group[i];
		if(vectorContainsElementAlready(&useSeqs, seq2) ==false)
		  {
		    dist = getDistance(seq1,seq2);
		    
		    /* using the similarity masking idea to reduce the number of potential seeds */
		    if(dist < distMaskCutoff)
		      {
			int len1 = alignPtr->getSeqLength(seq1 + 1); 
			int len2 = alignPtr->getSeqLength(seq2 + 1);
			
			if(len1 >= len2){useSeqs.push_back(seq2);}
			if(len1 <  len2){useSeqs.push_back(seq1);}
		      }		  
		    squareDist += pow(dist,2);	
		  }
	      }
	    
	    dist = sqrt(squareDist);
	    
	    seqList.push_back(seq1);
	    seqDist.push_back(dist);
	    
	  }
	
	shellSort(&seqList, &seqDist);
	
	
	
	int furthestObject  = seqList[seqList.size() - 1];
	
	seqDist.clear();     
	seqList.clear();
	
	return furthestObject;    
      }
    
    catch(const exception &ex)
      {
	cerr << "Exception in findFurthestObjectFromGroup function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }        
  } 
  void SparseMap::addSequenceAsSeed(int sequence, double distMaskCutoff)
  {
    try
      {
	bool addSequence = true;
	for(int i = 0; i < seedSeqs.size(); i++)
	  {
	    double dist = getDistance(sequence, seedSeqs[i]);
	    if(dist <= distMaskCutoff){addSequence = false;}
	  }
	if(addSequence == true){ seedSeqs.push_back(sequence);}
      }
    
    catch(const exception &ex)
      {
	cerr << "Exception in addSequenceAsSeed function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }    
  }
  void SparseMap::initialiseCoordinates()
  {
    try
      {
        data->clear();
        data->resize(numSeqs, Point(dimensions));
      }
    catch(const exception &ex)
      {
        cerr << "Exception in initialiseCoordinates function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }  
  void SparseMap::generateCoordinates(int numOutputDims, bool useSeedsOnly)
  {
    try
      {
	if(numOutputDims >= 1)
	  {
	    int numInputDims = dimensionGroups.size();
	    cout << "currently have " << numInputDims <<  " ref groups ";
	    cout << "but we WANT only the first " << numOutputDims << "\n";
	    
	    for(int dim = numInputDims; dim > numOutputDims; dim--)
	      {
		dimensionGroups.pop_back();
	      }
	    
	  }
	
	double dist;
	dimensions = dimensionGroups.size();
	
        //double sqrtDims = sqrt(dimensions);
	
	cout << "embedding " << numSeqs << " sequences into "<< dimensions <<" dimensions...\n";
	initialiseCoordinates();
	clock_t timeStart, timeStop; // cpu timer just for the embedding into dimensionGroups.
	time_t  timeBegin, timeEnd;  // wall timer for the same
	
	timeStart = clock();
	timeBegin = time(NULL);
	
	
	
	if(useSeedsOnly == true)
	  {
	    /* a number of advantages become apparent when seed sequences are used as reference groups. These are described
	       below (called FIXes). This is why we have a separate loop
	    */
	    
	    seedSeqs.clear(); 
	    for(int dim = 0; dim < dimensions; dim++)
	      {
		int seed = dimensionGroups[dim][0];			
		seedSeqs.push_back(seed);
		
	      }
	    // ok, now we have a list (again) of the seed sequences used.
	    
	    
	    /*
	      FIX1: we're switching the loop!!

	      The reason for this is simple. By switching the loop to
	      embed each sequence in turn into all dimensions, we make
	      no difference to the number of calculations. However, it
	      means that a full array does not need to be computed at
	      the one time. By simply storing the seeds, we can embed
	      virtually unlimited numbers of sequences, writing out
	      each vector to a file, one at a time. If all sequences
	      can't be held in memory, that's ok. All you really need
	      is one sequence, and a list of seeds, to create an
	      embedding of that one sequence.
	      
	      It will mean, however, that the shortcut about
	      re-calculating distances (FIX2, below) will become
	      ineffective; it's a trade off between speed and
	      capacity. Nevertheless, the method should still scale
	      linearly with number of seeds.
	   
	      

	    */
	    

	    /*USE THIS LOOP FOR FIX 1*/
	    /*
	    for(int seq1 = 0; seq1 < numSeqs; seq1++)
	      {
	      //cout << "embedding sequence " << seq1 + 1 << " / " << numSeqs << " into " << dimensions << " dimensions";
	      // seq1 is the sequence in the sequence list.
		for(int dim1 = 0; dim1 < dimensions; dim1++)
		  {
		  /// dim1 is the dimension into which the sequence is being embedded
		  /// seq2 is the seed sequence corresponding to that dimension. 
		    
		    int seq2 = seedSeqs[dim1];	    
	    */

	    
	    /* CURRENT LOOP IS TO ALLOW FIX 2 */
	    
	    
	    for(int dim1 = 0; dim1 < dimensions; dim1++)
	      
	      {
		//cout << "embedding sequence " << seq1 + 1 << " / " << numSeqs << " into " << dimensions << " dimensions";
		

		
		for(int seq1 = 0; seq1 < numSeqs; seq1++)


		  {
		  /// dim1 is the dimension into which the sequence is being embedded
		  /// seq2 is the seed sequence corresponding to that dimension. 
		    
		    int seq2 = seedSeqs[dim1];
		    
		    
		    /*
		      FIX2: we're cutting down on distance calculations!!! 

		      If the sequence to be embedded is also a seed
		      sequence, then the distance between that
		      sequence and another seed might already be
		      available elsewhere in the array.  This is only
		      effective if we have at least a partial array in
		      memory. This was impossible when reference
		      groups contained more than one sequence.
		      
		      When small numbers of seeds are used, this means
		      little improvement in time. However, if 10,000
		      seeds are used, well, then this fix means that
		      something like 50,000,000 distance calculations
		      do not need to be done.
		    */
		    


		  
                    /// Has the distance between these two sequences been calculated already?
		    
		    if((*data)[seq1][dim1] == 0)
		      {
			/// if not, then is seq1 also a seed sequence? 

			if(vectorContainsElementAlready(&seedSeqs,seq1) == false)
			  {
			    /// if not then just do the distance calculation and embed the sequence

                            (*data)[seq1][dim1] = getDistance(seq1,seq2);
			  }
			

			else
			  {
			    
			    /// if so, then the distance might already be available elsewhere in the array.
			    /// First, need the location of seq1 in seedSeqs. This is called dim2
			    
			    int dim2 = vectorContainsElementAt(&seedSeqs,seq1);
			    

			    /// has the corresponding cell already been filled?
			  if((*data)[seq2][dim2] != 0)
			    {
			      (*data)[seq1][dim1] = (*data)[seq2][dim2];
			    }
			  else
			    {
                              (*data)[seq1][dim1] = getDistance(seq1,seq2);
			    }
			  }
		      }
		  }
		//timeStop = clock();
		//cout  << " ... "<< (double)(timeStop-timeStart)/CLOCKS_PER_SEC << " s\n";
	      }
	    
	  }
	
	else if(useSeedsOnly == false)
	  {
	    for(int dim = 0; dim < dimensions; dim++)
	      {
		cout << "dim ("<< dim + 1 <<"/"<< dimensions<<") ("<<dimensionGroups[dim].size() <<" elements)...";
		
		for(int seq = 0; seq < numSeqs; seq++)
		  {
                    (*data)[seq][dim] = getMinDistanceInDimension(seq, dim);
		  }
		
		//timeStop = clock();
		//cout  << " ... "<< (double)(timeStop-timeStart)/CLOCKS_PER_SEC << " s\n";
		
	      }
	  }
	timeStop = clock();
	timeEnd = time(NULL);
	
	seedSeqs.clear();
	dimensionGroups.clear();
	useSeqs.clear();
	unsigned long int totalPossible = (((numSeqs * numSeqs) - numSeqs) / 2);

	cout << "\nSparseMap cpuEmbedTime:  "  << (double)(timeStop-timeStart)/CLOCKS_PER_SEC;
	cout << "\nSparseMap wallEmbedTime: "  << timeEnd - timeBegin;
	cout << "\nSparseMap Dimensions: "     << dimensions;
	cout << "\nSparseMap DistCalls:  "     << counter;
	cout << "\nFullMatrix DistCalls: "     << totalPossible;
	cout << "\n";
      }
    catch(const exception &ex)
      {
	cerr << "Exception in generateCoordinates function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }    
  }
  double SparseMap::getMinDistanceInDimension(int seq, int dim)
  {
    // This is a SparseMap function, not SeedMap. set min to be the first distance. 
    try
      {
        double min;
        // need to find a minimal number of dimensions for which to calculate the true minimal distance
        // in order to become comparable to the full distance matrix...This may take some time, but should
        // start on the test cases that display greatest difference between the Full and Embedded data,
        // and try to shorten the gap. 
	
        if(dim == 0) // if it is the first dimension, calculate all
	  {
            min = getDistance(seq, dimensionGroups[dim][0]);
            double dist;
            for(int element = 1; element < dimensionGroups[dim].size(); element++)
	      {
                dist = getDistance(seq, dimensionGroups[dim][element]);
                if(dist < min)
		  {
                    min = dist;
		  }
	      }        
	  }
        else
	  {
            // find the shortest euclidean distance(in dim-1 dimensions) 
            // and then get the real distance.
            int seqWithShortestApprox;
            seqWithShortestApprox = getSeqWithShortestApproxDistance(seq, dim);
            min = getDistance(seq, seqWithShortestApprox);
	  }
	
        return min;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in getMinDistanceInDimension function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }
  int SparseMap::getSeqWithShortestApproxDistance(int seq, int dim)
  {
    // This is a SparseMap function, not SeedMap. 
    // get euclidean distance between this one and all the others.
    try
      {
        double smallestDistSoFar;
        int indexOfSmallest;
        smallestDistSoFar = getApproxEuclideanDist(seq, dimensionGroups[dim][0], dim);
        indexOfSmallest = dimensionGroups[dim][0];
	
        double dist;
        for(int element = 1; element < dimensionGroups[dim].size(); element++)
	  {
            dist = getApproxEuclideanDist(seq, dimensionGroups[dim][element], dim);
            if(dist < smallestDistSoFar)
	      {
                indexOfSmallest = dimensionGroups[dim][element];
	      }
	  }
        return indexOfSmallest;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in getSeqWithShortestApproxDistance function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }
  double SparseMap::getApproxEuclideanDist(int seq1, int seq2, int dim) 
  {
    /* This is a SparseMap function, not SeedMap. 
       This is the Euclidean distance function to find the minimal distance
       between a point and a reference group. To slightly reduce the overhead
       of this function, use the squared Euclidean distances instead */
    try
      {
	double sumOfSquares = 0.0;
	for(int i = 0; i < dim; i++)
	  {
	    sumOfSquares += square((*data)[seq1][i] - (*data)[seq2][i]);
	  }
	
	return sumOfSquares;
      }
    catch(const exception &ex)
      {
	cerr << "Exception in getApproxEuclideanDist function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }    
  }
  void SparseMap::calculateDistanceMatrix()
  {
    try
      {
	
        cout << "calculating distance matrix\n";
        for(int i = 0; i < dimensionGroups.size(); i++)
	  {
            dimensionGroups[i].clear();
	  }
        dimensionGroups.clear();
	
        double dist;
        sparseMapDistMat.reset(new SymMatrix);
        sparseMapDistMat->ResizeRect(numSeqs + 1);
        
        for(int i = 0; i < numSeqs; i++)
	  {
            for(int j = i + 1; j < numSeqs; j++)
	      {
                dist = getEmbeddedDistance(&(*data)[i], &(*data)[j]);
                (*sparseMapDistMat)(i + 1, j + 1) = dist;
	      }
	  }
      }
    catch(const exception &ex)
    {
      cerr << "Exception in calculateDistanceMatrix function\n";
      cerr << ex.what() << endl;
      cerr << "Terminating program. Cannot continue\n";
      exit(1);
    }          
  }
  double SparseMap::getEmbeddedDistance(Point* vec1, Point* vec2)
  {
    
    /*
      this Euclidean distance function is used to create the
      distance matrix from the full set of coordinates. This
      function is defined separately from the getApproxEuclideanDist
      function (above) because of the different data structures used.

      NOTE: The distance function is initially the normal Euclidean measure,
      but to derive the distance from embedded coordinates, this value is 
      divided by the square root of the dimensionality of the vectors. 

    */
    try
      {
	
	
        double dist;
        double sumSquares = 0.0;
        int vec1Size = vec1->getNumDimensions();
        int vec2Size = vec2->getNumDimensions();
	
        int p = 2;       // p describes the Minkowski distance formula for Euclidean space
        int k = vec1Size;// k describes the dimensionality of the space.
	
        // when p == 2, then the p-root of a value  is equal to sqrt(value). Use this until we choose non-Euclidean Space!!
	
	
        for(int i = 0; i < vec1Size && i < vec2Size; i++)
          {
            sumSquares += square( ((*vec1)[i] - (*vec2)[i]) / sqrt(k) )    ;
          }
        dist = sqrt(sumSquares);
        return dist;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in getEmbeddedDistance function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }
  }
  double SparseMap::getEmbeddedDistance2(Point* vec1, Point* vec2, int k)
  {
    /*  
	this Euclidean distance function is used to create the
	distance matrix from a partially completed set of coordinates. This
	function is defined separately from the getApproxEuclideanDist
	function (above) because of the different data structures used.
    */
    try
      {


        double dist;
        double sumSquares = 0.0;
        int vec1Size = vec1->getNumDimensions();
        int vec2Size = vec2->getNumDimensions();

        int p = 2;       // p describes the Minkowski distance formula for Euclidean space
	
	
        // when p == 2, then the p-root of a value  is equal to sqrt(value). Use this until we choose non-Euclidean Space!!
	
	
        for(int i = 0; i < k; i++)
          {
            sumSquares += square( ((*vec1)[i] - (*vec2)[i]) / sqrt(k) )    ;
          }
        dist = sqrt(sumSquares);
        return dist;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in getEmbeddedDistance function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }
  }
  void SparseMap::printDistanceMatrixToFile()
  {
    //cout << "*****************************************\n";
    //cout << "** Printing dist matrix to distMat.out **\n";
    //cout << "*****************************************\n";
    
    ofstream outfile("distMat.out", ofstream::trunc);
    cout <<"writing distance matrix to distMat.out...";
    for(int i = 0; i < numSeqs; i++)
      {
        outfile << alignPtr->getName(i + 1) << "\t";
        for(int j = 0; j <= i; j++)
	  {
            outfile << setprecision(10) << (*sparseMapDistMat)(i + 1, j + 1) << " ";
	  }
        outfile << "\n";
      }  
    cout << "done\n";  
  }
  void SparseMap::printCoordinatesToFile()
  {
    //cout << "*********************************************\n";
    //cout << "** Printing coordinates to coordinates.out **\n";
    //cout << "*********************************************\n";
    ofstream outfile("coordinates.out", ofstream::trunc);
    cout << "writing vector data to coordinates.out...";    
    for(int i = 0; i < data->size(); i++)
      {
        outfile << alignPtr->getName(i + 1);
        for(int j = 0; j < (*data)[i].getNumDimensions(); j++)
	  {
            outfile  << "\t" << setprecision(10) <<  (*data)[i][j];
	  }
        outfile << "\t$\n";
      }
    cout << "done\n";
  }
  double SparseMap::getDistance(int seq1, int seq2) {

    try {
      dist = 0.0;
      counter2++;
      if(seq1 == seq2){
	dist = 0.0;
	
	counter4++;       
      }
      
      else { /* sequences not same label */
	
	/* Hash
	   generate hash,
	   first, exploit symmetry and order labels min/max
	   second, turn label pair into one number index=max*N+min
	   third, use golden ratio to map index on hash bucket
	*/
	int iMax, iMin;
	if (seq1 > seq2){
	  iMax = seq1+1;
	  iMin = seq2+1;
	}
	else{ 
	  iMax = seq2+1;
	  iMin = seq1+1;
	}
	int iIndex = numSeqs*iMax + iMin;
	iIndex = (iIndex > 0) ? iIndex : -iIndex;
#define GOLDEN_RATIO 0.6180339887498948482
	double dIndex = GOLDEN_RATIO * iIndex;
	dIndex  = fmod(dIndex, 1.00);
	dIndex *= iHashKey;
	iIndex  = (int)(dIndex);

	int iCount = 0;
	dist = -1.00;
	
	if (NULL == pprHash[iIndex]){
	  iCount = 0;
	}
	else {
	  for (iCount = 0; 1 != pprHash[iIndex][iCount].last; 
	       iCount++){
	    
	    if ( (pprHash[iIndex][iCount].first == iMax) && 
		 (pprHash[iIndex][iCount].secnd == iMin) ){
	      
	      dist = pprHash[iIndex][iCount].dist;
	      
	    } /* distance already calculated */
	    
	  } /* went through hash bucket, but not last element */
	  
	  /* have to do last element as well */
	  if ((pprHash[iIndex][iCount].first == iMax) &&
	      (pprHash[iIndex][iCount].secnd == iMin) ){
	    
	    dist = pprHash[iIndex][iCount].dist;
	    
	  } /* (last) distance already calculated */
	  else {
	    iCount++;
	  }
	  
	} /* hash bucket existed already */
	
	if (dist < 0.00){
	  /* at this stage we either have encountered the pair 
	     for the first time or there is a conflict */
	  
#ifdef DIAGNOSTIC
            if (iCount < MAX_COLLISIONS){
              ziCollisions[iCount]++;
            }
#endif
	  if (NULL == (pprHash[iIndex] = (myHash *)
		       realloc(pprHash[iIndex], 
			       (iCount+1)*sizeof(myHash)))){
	    printf("%s:%d: couldn't realloc hash bucket %d to %d\n", 
		   __FILE__, __LINE__, iIndex, iCount+1);
	    exit(-1);
	  }

	  if(fast) {
	    //FastPairwiseAlign fastAlign;
	    /* NOTE: 
	       this isn't the recommended function to use, however ... */
	    fastAlign.pairwiseAlign(alignPtr, &dist, seq1+1, seq2+1);
	    
	  } /* fast */
	  else {
	    //FullPairwiseAlign fullAlign;
	    
	    /* NOTE: 
	       this is not the recommended function to use, however ... */
	    fullAlign.pairwiseAlign(alignPtr, &dist, seq1+1, seq2+1);
	  } /* full */
	  
	  pprHash[iIndex][iCount].first = iMax; 
	  pprHash[iIndex][iCount].secnd = iMin;
	  pprHash[iIndex][iCount].dist  = dist;
	  pprHash[iIndex][iCount].last  = 1;
	  if (iCount > 0){
	    pprHash[iIndex][iCount-1].last  = 0;
	  }
	} /* had to calculate new distance */
	  
	counter++;
	sumDist+=dist;
	
      } /* sequences had different label */
      return dist;
    }
    catch(const exception &ex) {
      cerr << "Exception in getDistance function\n";
      cerr << ex.what() << endl;
      cerr << "Terminating program. Cannot continue\n";
      exit(1);
    }    

  } /* end of getDistance() */  



  int SparseMap::findLongestSequence()
  {    
    /* quite simply, find the longest sequence in the input sequences */
    try
      {                     
        int longestLen = 0;   
	int longestSeq = 0;
	
	for(int Seq=0; Seq < numSeqs; Seq++)
	  {
	    int Len = (int)(alignPtr->getSeqLength(Seq + 1));
	    if(Len > longestLen)
	      {
		longestLen = Len; longestSeq = Seq;	
	      }
	  }
	
        return longestSeq;    
      }
    
    catch(const exception &ex)
      {
        cerr << "Exception in  findLongestSequence function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }        
  }
  bool SparseMap::vectorContainsElementAlready(vector<int>* group, int element)
  {
    /* return true if a vector contains an element. 
       By default, set to work with ints...may need another if I'm working with doubles
    */
    try
      {    
        for(int i = 0; i < group->size(); i++)
	  {
            if((*group)[i] == element)
	      {
                return true;
	      }
	  }
        return false;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in vectorContainsElementAlready function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }
  int  SparseMap::vectorContainsElementAt(vector<int>* group, int element)
  {
    /* returns the index of the vector position of a particular element */
    try
      {    
        int j = 0;
	for(int i = 0; i < group->size(); i++)
	  {
            if((*group)[i] == element)
	      {
                j = i;
	      }
	  }
        return j;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in vectorContainsElementAt function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }    
  }
  void SparseMap::shellSort(vector<int>* seqList, vector<double>* seqVals)
  {
    /* 
       this is the start of the shellSort algorithm (that I adapted for two vectors...)
       This shellSort takes two vectors;
       
       seqList is a vector containing indices of the sequences themselves. They appear as they do in the input file
       seqVals is a vector containing values by which you want the sequences sorted (in ascending order)
    */
    try
      {
	
        int i, j, increment;
        double temp1;    // this is the value
        int temp2;    // this is the sequence itself 
        int numVals = seqList->size();
	
        increment = 3;
        while(increment > 0)
	  {
            for(i=0; i < numVals; i++)
	      {
                j = i;
                temp1 = (*seqVals)[i];
                temp2 = (*seqList)[i];
		
                while((j >= increment) && ((*seqVals)[j-increment] > temp1))
		  { 
                    (*seqVals)[j] = (*seqVals)[j - increment];
                    (*seqList)[j] = (*seqList)[j - increment];
                    j = j - increment;
		  }
		(*seqVals)[j] = temp1;     
                (*seqList)[j] = temp2;
	      }
            if(increment/2 != 0)
	      {
                increment = increment/2;
	      }
            else if(increment == 1)
	      {
                increment = 0;
	      }
            else
	      {
		increment = 1;
	      }
	  }
	
	
	
      }
    catch(const exception &ex)
      {
	cerr << "Exception in shellSort function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }
  }
  void SparseMap::evaluateEmbedding()
  {
    try
      {
	cout << "Evaluating embedding against full matrix in terms of\n";
	cout << "'goodness' & 'badness' of fit, using correlation & stress\n";	
	vector<double> X; // list to hold the original distances...
	vector<double> Y; // list to hold the embedded distances...


	for(int i = 0; i < numSeqs; i++)
	  {
	    for(int j = i + 1; j < numSeqs; j++)
	      {
                double x = getDistance(i,j);                              // get the real dis-similarity value
                double y = getEmbeddedDistance(&(*data)[i], &(*data)[j]); // the embedded distance value

                X.push_back(x);
                Y.push_back(y);
	      }
	  }
	cout << "\n";
	
	double correlation = calculateCorrelation(X,Y);
	double stress      = calculateStress(X,Y);

        
	cout << "SparseMap Stress:\t"      << stress      <<"\n";
	cout << "SparseMap Correlation:\t" << correlation << "\n";
	
	cout << "\n";
	
	X.clear();
	Y.clear();
      }
    catch(const exception &ex)
      {
	cerr << "Exception in evaluateEmbedding function\n";
	cerr << ex.what() << endl;
	cerr << "Terminating program. Cannot continue\n";
	exit(1);
      }
  }
  void SparseMap::evaluateEmbedding2()
  {
    try
      {
        cout << "Evaluating embedding against full matrix in terms of\n";
        cout << "'goodness' & 'badness' of fit, using correlation & stress\n";
        cout << "Evaluation is done across increasing values for k\n";
        cout << "Results are written to evalMat.out'\n";

        ofstream outfile("evalMat.out", ofstream::trunc);
        outfile << "Dims\tStress\tCorrelation\n";

        vector<double> X; // list to hold the original distances...
        vector<double> Y; // list to hold the embedded distances...

        int numDims = (*data)[0].getNumDimensions();

        X.clear();
        for(int i = 0; i < numSeqs; i++)
          {
            for(int j = i + 1; j < numSeqs; j++)
              {
                double x = getDistance(i,j); // get the real dis-similarity value

                X.push_back(x);
              }
          }

        for(int k = 1; k <= numDims; k++)
        {
            Y.clear();
            for(int i = 0; i < numSeqs; i++)
            {
                for(int j = i + 1; j < numSeqs; j++)
                {
                    double y = getEmbeddedDistance2(&(*data)[i], &(*data)[j],k); // the embedded distance value
                    Y.push_back(y);
                }
            }

            double correlation = calculateCorrelation(X,Y);
            double stress      = calculateStress(X,Y);
	    
            outfile << k << "\t" << setprecision(5) << stress << "\t" << setprecision(5) << correlation <<"\n";
        }

        X.clear();
        Y.clear();
      }
    catch(const exception &ex)
      {
        cerr << "Exception in evaluateEmbedding2 function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      }
  }
  double SparseMap::calculateCorrelation(vector<double> X, vector<double> Y)  
  {
    try
      {
	double xBar = 0.0;
	double yBar = 0.0;
	double xStd = 0.0;
	double yStd = 0.0;
        double correlation;
	
	int numElements = X.size();
	
	for(int i = 0; i < numElements; i++)
	  {
	    xBar += X[i];
	    yBar += Y[i];
	  }

	xBar = xBar / numElements;
	yBar = yBar / numElements;
	
	

	for(int i = 0; i < numElements; i++)
	  {
	    xStd += pow((X[i] - xBar),2);
	    yStd += pow((Y[i] - yBar),2);
	  }

	xStd = sqrt(xStd / (numElements - 1));
	yStd = sqrt(yStd / (numElements - 1));
	
	
        correlation = 0.0;
	for(int i = 0; i < numElements; i++)
	  {
	    correlation += (X[i]-xBar)*(Y[i]-yBar);
	  }
	
	correlation = correlation / ((numElements - 1) * xStd * yStd);

	return correlation;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in calculateCorrelation function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      } 
  }
  double SparseMap::calculateStress(vector<double> X, vector<double> Y)
  {
    try
      {
        double rStress  = 0.0; // the raw stress value; pretty uninformative.
        double nStress  = 0.0; // the normed stress value
        double kStress  = 0.0; // the Kruskal stress value

	int numElements = X.size();

        double  sum_dXY = 0.0;
	
	for(int i = 0; i < numElements; i++)
	  {
            rStress += pow((Y[i]  - X[i]),2);
            sum_dXY    += pow(X[i], 2);
	  }
	
        nStress = rStress / sum_dXY;

        kStress = sqrt(nStress);
    
        return kStress;
      }
    catch(const exception &ex)
      {
        cerr << "Exception in calculateStress function\n";
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
      } 
  }  
//  double SparseMap::getEuclideanDist(Point* vec1, Point* vec2)
//  {
//  /*
//      (IS THIS FUNCTION STILL BEING USED ANYWHERE???)
//
//      this Euclidean distance function is used to create the
//      distance matrix from the full set of coordinates. This
//      function is defined separately from the getApproxEuclideanDist
//      function (above) because of the different data structures used.
//  */
//    try
//      {
//
//        double dist;
//        double sumSquares = 0.0;
//        int vec1Size = vec1->getNumDimensions();
//        int vec2Size = vec2->getNumDimensions();
//        for(int i = 0; i < vec1Size && i < vec2Size; i++)
//	  {
//            sumSquares += square((*vec1)[i] - (*vec2)[i]);
//	  }
//        dist = sqrt(sumSquares);
//        return dist;
//      }
//    catch(const exception &ex)
//      {
//        cerr << "Exception in getEuclideanDist function\n";
//        cerr << ex.what() << endl;
//        cerr << "Terminating program. Cannot continue\n";
//        exit(1);
//      }
//  }

}
