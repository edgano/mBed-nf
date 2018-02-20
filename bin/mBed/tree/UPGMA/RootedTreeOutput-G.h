public:
	void printPhylipTreeNoDists(RootedGuideTree* tree, ofstream* ptrToFile,
                                    Alignment *alignPtr); /* FS, 2009-04-22 */

private:
	/* FS, 2009-04-22 --> */
	void phylipTraverseNoHeights(ofstream* ptrToFile, Alignment *alignPtr, Node* t);
	bool isLeafNode(Node* t); /* <-- FS */ 
