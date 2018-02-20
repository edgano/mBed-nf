

public:
	/* FS, 2009-04-22 --> */
	auto_ptr<AlignmentSteps> generateTree(RootedGuideTree* phyTree,
					      DistMatrix* distMat, SeqInfo* seqInfo,
					      bool overwrite, vector<int>* realNums,
					      ofstream* tree = 0); /* <-- FS */

private:
	/* FS, 2009-04-22 --> */
	Node **initialiseNodes(double *distanceMatrix, int firstSeq, vector<int>* realNums,
                               int nSeqs); /* <-- FS */
