
public:
	/* FS, 2009-04-22 --> */
	auto_ptr<AlignmentSteps> 
	  treeFromDistMatrix(RootedGuideTree* phyloTree,
			     DistMatrix* distMat, Alignment *alignPtr, int seq1,
			     int nSeqs, string& phylipName,
			     vector<int>* realNums, bool printTree = true); /* <-- FS */
