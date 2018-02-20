/**
 * Author: Gordon Blackshields
 *
 * for comments contact fabian.sievers@ucd.ie (FS)
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream>
#include <ctime>
#include "alignment/Alignment.h"
#include "alignment/Sequence.h"
#include "general/clustalw.h"
#include "general/UserParameters.h"
#include "substitutionMatrix/SubMatrix.h"
#include "general/Utility.h"
#include "fileInput/FileReader.h"
#include "interface/InteractiveMenu.h"
#include "interface/CommandLineParser.h"
#include "general/DebugLog.h"
#include "general/ClustalWResources.h"
#include "general/Stats.h"
#include "SparseMap.h" /* G */
#include "SmartKMeans.h" /* G */


namespace clustalw
{ 
  UserParameters* userParameters;	/* C & G */
  Utility* utilityObject;		/* C */
  SubMatrix *subMatrix; 		/* C & G */
  DebugLog* logObject;			/* C */
  Stats* statsObject;			/* C */
}
using namespace std;
using namespace clustalw;

void help(char *argv0) /* included argc0, FS, 2009-08-06 */
  // Help message for helpless users
{
  cout << "\n\ntypical usage: " << argv0 << " [-infile pathToInputFile] ";
  //cout << "[-seedfile pathToSeedFile] ";
  cout << "[-speed speed] ";
  cout << "[-clustering clusteringType] ";
  cout << "[-method methodType] ";
  cout << "[-numInputSeeds number] ";
  cout << "[-numOutputDims number] ";
  cout << "[-useSeedsOnly true] ";
  cout << "[-findOutliers true] ";
  cout << "[-seedPosition vector] ";
  cout << "[-seedSelection method] ";
  cout << "[-seedTreatment method] ";
  cout << "[-evaluate False]";
  cout << "[-evaluate2 False]";

  cout << " \n\nwhere:\n";
  cout << " [-speed] can be either fast or slow (default = fast)\n";
  cout << " [-clustering] can be either upgma, kmeans, hybrid or none (default = upgma).\n";
  cout << " [-method] can be either SparseMap or SeedMap or FullMatrix\n (default = FullMatrix).\n";
  cout << " [-numInputSeeds] is the number of initial seed sequences to use to start the pre-processing (SeedMap).\n";
  cout << " [-numOutputDims] is the final number of output dimensions (SeedMap, SparseMap).\n";
  cout << "                  It embeds the sequences with respect to the first k reference groups. \n";
  cout << " [-useSeedsOnly] is 'true' by default; forces the program to embed all sequences using just the seeds.\n";
  cout << " [-findOutliers] is 'true' by default; setting it as 'false' will gain speed (but may lose accuracy.)\n";
  cout << " [-seedfile] allows the user to supply a set of sequences (in a FASTA format file) to use as seeds.\n";
  cout << " [-seedSelection] may be 'naive','single', 'medoids' or 'multi'. This defines how the seeds are chosen.\n";
  cout << " [-seedSelection] may be 'naive','bySeqLen','bySeqDist' or 'bySeqFeats'. This defines how the seeds are treated.\n";
  cout << " [-evaluate] allows you to test the resultant vectors by correlation with the orginal distances\n";
  cout << " \n\nIf you do not specify any options, then the program will create a full distance matrix and upgma tree by default\n";
  exit(1);
}



int main(int argc, char *argv[]){      

  userParameters = new UserParameters(false);	/* C & G */
  utilityObject = new Utility();		/* C */
  subMatrix = new SubMatrix();			/* C */
  statsObject = new Stats();			/* C */
  ClustalWResources *resources = ClustalWResources::Instance(); /* C */
  resources->setPathToExecutable(string(argv[0])); /* C */
  userParameters->setDisplayInfo(true); /* C */
  vector<string> args; /* FS, 2009-04-21 */
  
  /* G --> */

  string inputFile;
  string seedFile;
  string speed = "fast";
  string clustering = "upgma";
  string method = "FullMatrix";

  string seedSelectionMethod = "uniform";
  string seedTreatmentMethod = "uniform";
  
  int numInputSeeds = -1;
  int numOutputDims = -1;
  double distMaskCutoff = 0.0;
  
  bool fast = true;
  bool sparse = false;
  bool hybrid = false;
  bool upgma = true;
  bool kmeans = false;
  bool doClustering = true;
  bool useSeedsOnly = true;
  bool findOutliers = true;
  bool useSeedFile = false;
  bool evaluate = false;
  bool evaluate2 = false;
  /* <-- G */
  
  
  //userParameters->setDebug(5);       
#if DEBUGFULL    
  if(DEBUGLOG) {
    cout << "debugging is on\n\n\n";
    logObject = new DebugLog("logfile.txt");
    logObject->logMsg("Loggin is on!");
  }
#endif

  cout << endl;
  cout << "+-----------------------------------------+" << endl;
  cout << "| mBed -- fast clustering / tree building |" << endl;
  cout << "+-----------------------------------------+" << endl;
  cout << "Author: Gordon Blackshields, UCD\n" << endl;
  cout << "for help type " << argv[0] << " -help\n" << endl;
  if (1 == argc){
    exit(1);
  }

  for (int i = 1; i < argc; ++i) {

    if (strcmp(argv[i], "-help") == 0) {
      help(argv[0]);
      exit(1);
    } /* help */
    
    else if(strcmp(argv[i], "-infile") == 0)
      {
	inputFile = string(argv[i+1]);
	i++; /* FS, 2009-04-21 */
	
	bool flag = false;
	ifstream fin;
	fin.open(inputFile.c_str(), ifstream::in);
	
	if( fin.is_open() ){
	  flag=true;
	}
	
	fin.close();
	
	if (flag == true) 
	  {
	    cout << "using input file " << inputFile <<"\n";
	  }
	else if (flag == false) 
	  {
	    cout << "error: input file does not exist: "<<inputFile <<"\n";
	    exit(1);
	  }
      } /* file name */
    
    else if (strcmp(argv[i], "-seedfile") == 0)
      {
	
	seedFile = string(argv[i+1]);
	bool flag = false;
	i++; /* FS, 2009-04-21 */
	
	ifstream fin;
	fin.open(seedFile.c_str(), ifstream::in);
	if( fin.is_open() ) 
	  {
	    flag=true;
	  }
	
	fin.close();
	if (flag == true) 
	  {
	    cout << "using seed file " << seedFile <<"\n";
	    useSeedFile = true;
	  }
	else if (flag == false) 
	  {
	    cout << "error: seed file does not exist: "<<seedFile <<"\n";
	    exit(1);
	  }
      } /* seedfile */
//**************************************************//
//**************************************************//
                //seedPosition

//**************************************************//
//**************************************************//
    else if (strcmp(argv[i], "-method") == 0)
      {
	method = argv[i+1];
	cout << "Method:\t" << method <<"\n";
	i++; /* FS, 2009-04-21 */

	
	if(method == "SparseMap") 
	  {
	    sparse = true;
	  }

	else if(method == "SeedMap") 
	  {
	    sparse = true;
	  }
	else if(method == "FullMatrix") 
	  {
	    sparse = false;
	  }
      else if ( (method != "SparseMap")  &&
		(method != "SeedMap")    &&
		(method != "FullMatrix") ) 
	{
	  cout << "\nInvalid argument for " << argv[i] <<" command: " << method <<".\n";
	  cout << "Please specify one of 'SparseMap','SeedMap','FullMatrix'\n";
	  help(argv[0]);
	  exit(1);
	}
	//cout << method+"\n";
      } /* method */
    

    else if(strcmp(argv[i], "-speed") == 0)
      {
	
	speed = argv[i+1];
	i++; /* FS, 2009-04-21 */
	
	cout << "Speed:\t" << speed <<"\n";
	
	if(speed == "fast") 
	  {
	  fast = true;
	  }
	else if(speed == "slow") 
	  {
	    fast = false;
	  }
      else if ( (speed != "fast") &&
		(speed != "slow") ) 
	{
	  cout << "\nInvalid argument for " << argv[i] <<":\t"<<speed<<".\n";
	  cout << "Please specify 'fast' or 'slow'\n";
	  help(argv[0]);
	  exit(1); // multiple exit points in parsing algorithm
	}
      } /* speed */
    

    else if(strcmp(argv[i], "-clustering") == 0)
      {
	
	clustering = argv[i+1];
	cout << "Clustering:\t" << clustering <<"\n";
	i++; /* FS, 2009-04-21 */
	
	if(clustering == "upgma") 
	  {
	    upgma  = true;
	    hybrid = false;
	    kmeans = false;
	    doClustering = true;
	  }
	else if(clustering == "kmeans") 
	  {
	    upgma  = false;
	    kmeans = true;
	    hybrid = false;
	    doClustering = true;
	  }
	else if(clustering == "hybrid") 
	  {
	    upgma  = false;
	    kmeans = true;
	    hybrid = true;
	    doClustering = true;
	  }
	else if(clustering == "none") 
	  {
	    upgma  = false;
	    kmeans = false;
	    hybrid = false;
	    doClustering = false;
	  }
	else if ( (clustering != "upgma")  &&
		  (clustering != "kmeans") &&
		  (clustering != "hybrid") &&
		  (clustering != "none") ) 
	  {
	    cout << "\nInvalid argument for " << argv[i] << ":\t" << clustering <<".\n";
	    cout << "Please specify one of 'upgma','hybrid','kmeans','none'\n";
	    help(argv[0]);
	    exit(1);
	  }
      } /* clustering */
    
    else if(strcmp(argv[i], "-numInputSeeds") == 0)
      {
	
	numInputSeeds = atoi(argv[i+1]);
	i++; /* FS, 2009-04-21 */
	
	cout << "numInputSeeds:\t" << numInputSeeds << "\n";
	
	if(numInputSeeds < 1) 
	  {
	    cout << "Error: please specify a positive integer for number of input seeds\n";
	    help(argv[0]);
	    exit(1);
	  }
      } /* number of seeds */
 

    else if(strcmp(argv[i], "-distMask") == 0)
      {
	distMaskCutoff = atof(argv[i+1]);
	i++; /* FS, 2009-04-21 */
	
	cout << "using distance masking cutoff of:\t" << distMaskCutoff << "\n";
	
	if(distMaskCutoff > 1.0 ) 
	  {
	    cout << "Error: please specify a distance mask cutoff of between 0 and 1.\n";
	    help(argv[0]);
	    exit(1);
	  }
      } /* distance mask cutoff value */

   

    else if(strcmp(argv[i], "-numOutputDims") == 0)
      {
	
	numOutputDims = atoi(argv[i+1]);
	cout << "numOutputDims:\t" << numOutputDims << "\n";
	i++; /* FS, 2009-04-21 */

      if(numOutputDims < 1) 
	{
	  cout << "Error: please specify a positive integer for number of output dimensions\n";
	  help(argv[0]);
	  exit(1);
	}
      } /* dimensions */
    
    else if(strcmp(argv[i], "-useSeedsOnly") == 0)
      {
	
	string option = argv[i+1];
	cout << "Using input seeds as single-element reference groups:\t" << option << "\n";
	i++; /* FS, 2009-04-21 */
	
	if(option == "true") 
	  {
	    useSeedsOnly = true;
	  }
	else if(option == "false") 
	  {
	    useSeedsOnly = false;
	  }

      else if(option != "true" && option !="false") 
	{
	  cout << "\nInvalid argument for " << argv[i] <<":\t " << option <<".\n";
	  cout << "Please specify 'true' or 'false'\n";
	  help(argv[0]);
	  exit(1); 
	}
      } /* seed only */
    
    else if (strcmp(argv[i], "-findOutliers") == 0)
      {
	
	string option = argv[i+1];
	i++; /* FS, 2009-04-21 */
	
	if(option == "true") 
	  {
	    findOutliers = true;
	  }
	else if(option == "false") 
	  {
	    findOutliers = false;
	  }
	else if(option != "true" && option !="false") 
	  {
	  cout << "\nInvalid argument for " << argv[i] <<":\t " << option <<".\n";
	  cout << "Please specify 'true' or 'false'\n";
	  help(argv[0]);
	  exit(1); 
	  }
      } /* outliers */
    
    else if(strcmp(argv[i], "-seedSelection") == 0)
    {
      seedSelectionMethod = argv[i+1];
      cout << "seed Selection Method:\t" << seedSelectionMethod <<"\n";
      i++; /* FS, 2009-04-21 */
      if(
         (seedSelectionMethod != "minimal") &&
		 (seedSelectionMethod != "uniform") &&
		 (seedSelectionMethod != "length")  &&
         (seedSelectionMethod != "random") &&
		 (seedSelectionMethod != "medoid")
		 )
	{
	  cout << "unrecognised option for " << argv[i] << " command: " << seedSelectionMethod << "\n";
	  cout << "Please choose one of the following; minimal, uniform, length, medoid or random.\n";
	  help(argv[0]);
	  exit(1);
	}
    } /* seed selection */
    


    else if(strcmp(argv[i], "-seedTreatment") == 0)
      {
	
	seedTreatmentMethod = argv[i+1];
	i++; /* FS, 2009-04-21 */
	
	cout << "seed Treatment Method:\t" << seedTreatmentMethod <<"\n";
	
	if( (seedTreatmentMethod != "naive") && 
	    (seedTreatmentMethod != "none") && 
	    (seedTreatmentMethod != "usePivotObjects") &&
	    (seedTreatmentMethod != "usePivotGroups") && 
	    (seedTreatmentMethod != "usePivotObjectsAndConverge") &&
	    (seedTreatmentMethod != "usePivotGroupsAndConverge")  && 
	    (seedTreatmentMethod != "useSeqFeatures") )
	  {
	    cout << "unrecognised option for " << argv[i] << " command: " << seedTreatmentMethod << "\n";
	    cout << "Please choose one of the following; naive, usePivotObjects,usePivotGroups, usePivotObjectsAndConverge, usePivotGroupsAndConverge, useSeqFeatures\n";
	    help(argv[0]);
	    exit(1);
	  }
      } /* seed treatment */
  
  
    else if(strcmp(argv[i], "-evaluate") == 0)
      {
	
	string option = argv[i+1];
	cout << "you have chosen to evaluate the embedding (correlation, stress)\n";
	i++; /* FS, 2009-04-21 */
	
	if(option == "true") 
	  {
	    evaluate = true;
	  }
	else if(option == "false") 
	  {
	    evaluate = false;
	  }
	else if(option != "true" && option !="false") 
	  {
	    cout << "\nInvalid argument for " << argv[i] <<":\t " << option <<".\n";
	    cout << "Please specify 'true' or 'false'\n";
	    help(argv[0]);
	    exit(1); 
	  }
      }
    else if(strcmp(argv[i], "-evaluate2") == 0)
      {

        string option = argv[i+1];
        cout << "you have chosen to evaluate the embedding (correlation, stress)\n";
        cout << "across increasingly big values of k \n";
        i++; /* FS, 2009-04-21 */

        if(option == "true")
          {
            evaluate2 = true;
          }
        else if(option == "false")
          {
            evaluate2 = false;
          }
        else if(option != "true" && option !="false")
          {
            cout << "\nInvalid argument for " << argv[i] <<":\t " << option <<".\n";
            cout << "Please specify 'true' or 'false'\n";
            help(argv[0]);
            exit(1);
          }
      }
    
    else { /* FS, 2009-04-21 
	      arguments that are not recognised by Gordon's `logic`
	      are copied into the args vector to be processed by Clustal */

      args.push_back(argv[i]);

    } /* not seed map but clustal flag */

  } /* i <= 1 < argc */
  
    

  /* FS, 2009-04-21 
     Gordon's `logic` has already identified the infile
     let Clustal know about it */
  string zcAux;
  zcAux = string("-infile=");
  zcAux.append(inputFile);
  args.push_back(zcAux); 

  /*string zcAux2;          it will be used to add ref as file -> need to parse in the command line
  zcAux2 = string("-seedfile=");
  zcAux2.append(seedFile);
  args.push_back(zcAux2);*/


  cout << "Other Flags:" << endl;
  for (int i = 0; i < args.size(); i++){
    cout << i << ": " << args[i] << endl;
  }

  CommandLineParser cmdLineParser(&args, false);

  /*
    TIMER FUNCTIONS

    The first set are CPU timers - more accurate than wall times,
    but only go to a certain limit and then count backwards 
    (and go into negative numbers...)
    these timers are suffixed by "Start" or "Stop".

    The second set are the equivalent WALL timers.
    They are useful for long execution times, but only give integer values
    These timers are suffixed by "Begin" or "End"
  */
  
  clock_t readStart, readStop;           // clock to read in sequence data
  clock_t embeddingStart, embeddingStop; // clock for vector calculation
  clock_t matrixStart, matrixStop;       // clock for matrix calculation
  clock_t clusterStart, clusterStop;     // clock for cluster calculation.
    
  time_t readBegin, readEnd;             
  time_t embeddingBegin, embeddingEnd;
  time_t matrixBegin, matrixEnd;
  time_t clusterBegin, clusterEnd;
  
  double cpuReadTime;    int wallReadTime;
  double cpuEmbedTime;   int wallEmbedTime;
  double cpuMatrixTime;  int wallMatrixTime;
  double cpuClusterTime; int wallClusterTime;


  /* Timer 1: Time taken to read sequence file */
  readStart = clock(); readBegin = time(NULL);

  Alignment alignObj;	/* G */
  DistMatrix fullDistMat; /* G */
  FileReader myFileReader; /* F */
  myFileReader.seqInput(&alignObj, false);  /* F */

  /*FileReader mySeedReader;            will be the seeds file
  mySeedReader.seqInput(&alignObj, false);*/

  int numSeqs = alignObj.getNumSeqs();

  readStop = clock(); readEnd = time(NULL);
  
  cpuReadTime  = (double)(readStop - readStart)/CLOCKS_PER_SEC;
  wallReadTime = readEnd - readBegin;
  cout << "\nfileReadTime (s) [Wall / CPU]: " << wallReadTime << " " << cpuReadTime << "\n"; 
  

  /* G --> */
  SparseMap myMap(&alignObj, fast); 
  
  auto_ptr<DistMatrix> distMat;
  DistMatrix* ptrToMat;
  
  int nseqs = myMap.getNumSeqs();

  if(sparse == true) {

    vector<int> seeds;
    
    if(useSeedFile == true){
      /*
	if the seed file is supplied, then we need to open it and
	extract the sequences.  Right now it's very simple and naive:
	I'm not using the sequences from this seed file, but I'm
	searching through the alignment object from the original
	sequences, and recording the indices of these seeeds. This is
	just a temporary measure, as it was easier to code
	up. However, it also means that the seed sequences are needed
	to be present in the input file. To move away from this, I
	need to change the code to be able to use any sequences as
	seeds, regardless if they're in the sequence file or not.
      */
      
      /* <<CONTINUE HERE>> */


//*************************************************//
//      Generating the seed position vector        //
//*************************************************//

      for (int i=0; i<numInputSeeds;i++){
        seeds.push_back(i);
      }
//*************************************************//
//*************************************************//
    } /* (useSeedFile == true) */
    
    
    /*Timer 2: Time taken to do the embedding */
    embeddingStart = clock(); embeddingBegin = time(NULL);
    
    myMap.doEmbedding(method, numInputSeeds, numOutputDims, useSeedsOnly,
		      findOutliers, seeds, seedSelectionMethod, 
		      seedTreatmentMethod, distMaskCutoff);
    
    embeddingStop = clock(); embeddingEnd = time(NULL);
    cpuEmbedTime   = (double)(embeddingStop - embeddingStart)/CLOCKS_PER_SEC;
    wallEmbedTime  = embeddingEnd - embeddingBegin;
    
    cout << "\nEmbedding Time (s) [Wall / CPU]: " << wallEmbedTime << " " << cpuEmbedTime << "\n"; 
    
    //cout << "Embedding Time (s) [Wall / CPU]: " << embeddingEnd - embeddingBegin << " " <<(double)(embeddingStop-embeddingStart)/CLOCKS_PER_SEC << "\n";
    
    
    if (evaluate == true) 
      {
	myMap.evaluateEmbedding();
      }
    if (evaluate2 == true)
      {
        myMap.evaluateEmbedding2();
      }
    
    myMap.printCoordinatesToFile();
    
    if (doClustering == true)
      {
	
	if (upgma == true)
	  {
	    
	    
	    /*The actual clustering is done later - this loop is only for the distance matrix creation*/
	    /*Timer 3: Time taken to create a distance matrix (only valid for UPGMA clustering */
	    
	    matrixStart = clock(); matrixBegin = time(NULL);
	    myMap.calculateDistanceMatrix();
	    matrixStop = clock(); matrixEnd = time(NULL);
	    
	    cpuMatrixTime  = (double)(matrixStop - matrixStart)/CLOCKS_PER_SEC;
	    wallMatrixTime = matrixEnd - matrixBegin;
	    
	    cout << "\nMatrix Time (s) [Wall / CPU]: " << wallMatrixTime << " " << cpuMatrixTime << "\n"; 
	    //cout << "Matrix Time (s) [Wall / CPU]: " <<  matrixEnd - matrixBegin << " " << (double)(matrixStop-matrixStart)/CLOCKS_PER_SEC << "\n";
	    
	    
	    myMap.printDistanceMatrixToFile();
	    distMat = myMap.getDistMat();
	    ptrToMat = distMat.get();
	    
	  } /* (upgma == true) */
	
      
	else if ((kmeans == true) || (hybrid == true)) 
	  {
	    
	    
	    /*
	      Cluster the data with kmeans. 
	      
	      numIterations is the max number of iterations k-means will
	      have to find the optimal clustering (may terminate early).
	      
	      numCenters is the number of clusters k-means will find
	      
	    */
	    
	    int numCenters = 300; //nseqs / 50; /* <<DANGER>> */
	    if (numCenters <= 4) {
	      numCenters = 4;
	    }
	    int numIterations = (int)(log2(numCenters));
	    
	    auto_ptr<vector<Point> > data = myMap.getData(); 
	    auto_ptr<vector<Point> > thisClustering;
	    
	    /*Timer 4: time taken to cluster the objects (through whatever route */
	    clusterStart = clock(); clusterBegin = time(NULL);
	    
	    SmartKMeans kmeans;
	    thisClustering = kmeans.initClustering(data.get(), numCenters);
	    auto_ptr<KMeansResult> results;
	    results = kmeans.doClustering(data.get(), thisClustering.get(), numIterations);
	    
	    clusterStop = clock(); clusterEnd = time(NULL);
	    
	    if (hybrid == false)
	      {
		/*
		  That's it. Kmeans clustering only. UPGMA clustering will not 
		  be performed on the kmeans output 
		*/
		
		cpuClusterTime  = (double)(clusterStop - clusterStart)/CLOCKS_PER_SEC;
		wallClusterTime = clusterEnd - clusterBegin;
		
		cout << "\nClustering Time (s) [Wall / CPU]: " << wallClusterTime << " " << cpuClusterTime << "\n"; 
		
		
		
		vector<vector<int> > clusters;
		clusters.resize(numCenters);
		vector<int>* Assign = results->clusterAssignment.get();
		
		ofstream outfile("kmeans.out", ofstream::trunc);
		for(int i = 0; i < nseqs; i++){
		  int ClusterAssignment = (*Assign)[i];
		  
		  outfile << (*Assign)[i] << " ";
		}
		
		// this next part is unnecessary - it's purely for printing to the screen!! 
		/*
		  for(int i = 0; i < clusters.size(); i++) {
		  cout << "CLUSTER " << i << ": ";
		  for(int j = 0; j < Assign->size(); j++) {
		  if((*Assign)[j] == i) {
		  cout << alignObj.getName(j + 1) << " ";
		  
		  clusters[i].push_back(j);
		  }
		  }
		  cout << "\n";
		  }
		*/
		clusters.clear();
	  
	      } /* (hybrid == false) */
	    
	    else 
	      {
		// "hybrid" means to take the kmeans output and cluster with UPGMA
		cout << "\nKMEANS Time (s) [Wall / CPU]: " << wallClusterTime << " " << cpuClusterTime << "\n";
		TreeInterface tree;
		tree.doKMeansTree(&alignObj, data, results);
			      
		cpuClusterTime  = (double)(clusterStop - clusterStart)/CLOCKS_PER_SEC;
		wallClusterTime = clusterEnd - clusterBegin;
		
		cout << "\nClustering Time (s) [Wall / CPU]: " << wallClusterTime << " " << cpuClusterTime << "\n"; 
		
	      } /* !(hybrid == false) */

      } /* ((kmeans == true) || (hybrid == true)) */

    } /* (doClustering == true) */

  } /* (sparse == true) */



    if (sparse == false)
      {
	matrixStart = clock(); matrixBegin = time(NULL);
	fullDistMat.ResizeRect(nseqs + 1);
	double dist;
	if (fast == true) 
	  {	
	    
	    clustalw::FastPairwiseAlign fastAlign;
	    for(int i = 0; i < nseqs; i++) 
	      {
		cout << "computing matrix row ("<< i + 1 << "/" << nseqs << ")\n";
		for(int j = 0; j <= i; j++) 
		  {
		    if(i == j) 
		      {
			fullDistMat(i + 1, j + 1) = 0.0;
		      }
		    else 
		      {
			fastAlign.pairwiseAlign(&alignObj, &dist, i+1, j+1);
			fullDistMat(i + 1, j + 1) = dist;
		      }
		  }
	      }
	  } /* (fast == true) */
	
	else 
	  {
	    // use the slow Needleman-Wunsch distance measure
	    clustalw::FullPairwiseAlign fullAlign;
	    
	    for(int i = 0; i < nseqs; i++)
	      {
		cout << "computing matrix row ("<< i + 1 << "/" << nseqs << ")\n";
		for(int j = 0; j <= i; j++) 
		  {
		    if(i == j) 
		      {
			fullDistMat(i + 1, j + 1) = 0.0;
		      }
		    else 
		      {
			fullAlign.pairwiseAlign(&alignObj, &dist, i+1, j+1);
			fullDistMat(i + 1, j + 1) = dist;
		      }
			}
	      }
	  } /* !(fast == true) */
	
	
      
      /*
       * Get a pointer to the distance matrix to be used for clustering
       */
      
      ptrToMat = &fullDistMat;
      matrixStop = clock(); matrixEnd = time(NULL);
      cpuMatrixTime  = (double)(matrixStop - matrixStart)/CLOCKS_PER_SEC;
      wallMatrixTime = matrixEnd - matrixBegin;
      
      cout << "\nMatrix Time (s) [Wall / CPU]: " << wallMatrixTime << " " << cpuMatrixTime << "\n"; 
      long unsigned int totalPossible = (((nseqs * nseqs) - nseqs) / 2);
      cout << "\nFullMatrix DistCalls:  " << totalPossible << "\n";
      ofstream outfile("distMat.out", ofstream::trunc);
      for(int i = 0; i < nseqs; i++) 
	  {
		outfile << alignObj.getName(i + 1) << " ";
		for(int j = 0; j <= i; j++) 
		{
			outfile << setprecision(10) << fullDistMat(i + 1, j + 1) << " ";
		}
		outfile << "\n";
      }

    } /* (sparse == false) */
  

  if (upgma == true) 
  {
    // now that we have a distance matrix, apply UPGMA to cluster
    clusterStart = clock(); clusterBegin = time(NULL);
    cout << "we're at the cluster start point \n";
    
    TreeInterface tree;
    bool success = false;
    //string treeName = "temptree";
    string treeName = inputFile + ".dnd";
    //cout << "initialised the tree object \n";
    tree.doTreeOnly(ptrToMat, &alignObj, 1, nseqs, &treeName, &success);
    //cout << "and completed the tree object \n";
    clusterStop = clock(); clusterEnd = time(NULL);
    cpuClusterTime  = (double)(clusterStop - clusterStart)/CLOCKS_PER_SEC;
    wallClusterTime = clusterEnd - clusterBegin;

	cout << "\nClustering Time (s) [Wall / CPU]: " << wallClusterTime << " " << cpuClusterTime << "\n"; 
    
  } /* (upgma == true) */



    delete userParameters;
    delete utilityObject;
    delete subMatrix;
    
    if(logObject)
    {
        delete logObject;
    }


    cout << "\n"
      "+---------------------------+\n"
      "| Program ran to completion |\n"
      "+---------------------------+\n" << endl; 

    return 0;
}

