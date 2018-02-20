/* 
 * File: TreeInterface-G.h
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d in the last statement 
 *       of the TreeInterface class of TreeInterface.h, 
 *       that is, inside the clustalw namespace{} and the class TreeInterface{}
 *
 *       clustalw namespace{
 *       class TreeInterface{
 *          public:
 *           ...
 *          private:
 *           ...
 *       #include "TreeInterface-G.h"
 *       };
 *       }
 * Date: 2009-06-17 (FS)
 */

public:
	/* FS, 2009-03-27 --> */
	auto_ptr<RootedGuideTree> 
	  doTreeOnly(DistMatrix* distMat, Alignment *alignPtr, 
		     int seq1, int nSeqs, string* phylipName, bool* success);
	/* <-- FS */
	/* FS, 2009-04-22 --> */
	auto_ptr<RootedGuideTree> 
	  doKMeansTree(Alignment* alignObj, auto_ptr<vector<Point> > data,
		       auto_ptr<KMeansResult> results);
	/* <-- FS */

private:
	/* FS, 2009-04-22 --> */
	auto_ptr<AlignmentSteps> 
	  generateTreeFromDistMatUPGMA(RootedGuideTree* guideTree,
				       DistMatrix* distMat, Alignment *alignPtr, int seq1, int nSeqs,
				       string* phylipName, bool* success, bool printTree,
				       vector<int>* realNums = 0); 
	void printTree(RootedGuideTree* tree, int fSeq, Alignment* align); /* <-- FS */

