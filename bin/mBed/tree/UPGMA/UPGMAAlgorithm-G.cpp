


/* FS, 2009-04-22 --> */
auto_ptr<AlignmentSteps> UPGMAAlgorithm::generateTree(RootedGuideTree* phyTree,
						      DistMatrix* distMat,
						      SeqInfo* seqInfo, bool overwrite, vector<int>* realNums,
						      ofstream* tree)
{
  if (tree == 0 || !tree->is_open())
    {
      verbose = false;
    }
  //verbose = true;
  if (verbose)
    {
      cout << "\n\n\t\t\tUPGMA Method\n"
	   << "\n\n This is a ROOTED tree\n"
	   << "\n Numbers in parentheses are branch lengths\n\n";
    }

  progSteps.reset(new AlignmentSteps);

  Node** clusters;
  Node* root;
  numSeqs = seqInfo->numSeqs;
  const int sizeDistMat = ((numSeqs + 1) * (numSeqs + 2)) / 2;

    double* elements = overwrite ?
      distMat->getDistMatrix(seqInfo->firstSeq, seqInfo->numSeqs) :
      (double *)memcpy(new double[sizeDistMat],
		       distMat->getDistMatrix(seqInfo->firstSeq, seqInfo->numSeqs),
		       sizeDistMat * sizeof(double));

    clusters = initialiseNodes(elements, seqInfo->firstSeq, realNums, numSeqs);
    root = doUPGMA(clusters, tree);

    phyTree->setRoot(root);
    delete [] clusters;

    if(!overwrite)
      {
        delete [] elements;
      }
    distMat->clearSubArray();

    return progSteps;
}


Node** UPGMAAlgorithm::initialiseNodes(double* distanceMatrix, int fSeq,
                                       vector<int>* realNums, int nSeqs)
{
  int firstSeq = fSeq;

  bool useRealIndices = false;
  if(realNums != 0)
    {
      if(realNums->size() == nSeqs)
        {
	  useRealIndices = true;
        }
    }

  Node** nodes = new Node*[numSeqs];
  Node** nodeIter = nodes;

  *nodes = new Node(firstSeq, 0, 0);

  if(useRealIndices)
    {
      (*nodes)->realIndex = (*realNums)[0] + 1;
    }

  distanceMatrix++;

  // Move to first position in distanceMatrix.
  for(int elementIndex = 1, e = numSeqs; elementIndex < e;
      distanceMatrix += ++elementIndex)
    {
      Node* newcluster = new Node(elementIndex + firstSeq,
				  distanceMatrix, elementIndex);

      if(useRealIndices)
        {
	  newcluster->realIndex = (*realNums)[elementIndex] + 1;
        }
      (*nodeIter++)->next = newcluster;
      *nodeIter = newcluster;
    }

  return nodes;
}


