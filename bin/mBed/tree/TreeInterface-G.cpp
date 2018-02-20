/* 
 * File: TreeInterface-G.cpp
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d after the last function
 *       in TreeInterface.cpp but inside the clustalw namespace{} 
 *
 *       clustalw namespace{
 *         firstFunction(){}
 *           ...
 *         lastFunction(){}
 *       #include "TreeInterface-G.cpp"
 *       }
 * Date: 2009-06-17 (FS)
 */



/* FS, 2009-03-27 --> */
auto_ptr<RootedGuideTree> 
TreeInterface::doTreeOnly(DistMatrix* distMat, Alignment *alignPtr, int seq1, int nSeqs,
			  string* phylipName, bool* success)
{

  auto_ptr<RootedGuideTree> guideTree;
  guideTree.reset(new RootedGuideTree);
  generateTreeFromDistMatUPGMA(guideTree.get(), distMat, alignPtr, seq1, nSeqs,
			       phylipName, success);
  //  generateTreeFromDistMatUPGMA(guideTree.get(), &*distMat, &*alignPtr, (int)(seq1), (int)(nSeqs),
  //		       &*phylipName, &*success);
  return guideTree;
}
/* <-- FS */
/* FS, 2009-04-22 --> */
/**
 * The following function is used to build up a rooted tree from the KMeans clustering.
 *
 */
auto_ptr<RootedGuideTree> 
TreeInterface::doKMeansTree(Alignment* alignObj,
			    auto_ptr<vector<Point> > data, auto_ptr<KMeansResult> results)
{
  int NOCHILD = 0;

  // First print out the number of centers; This will be the size of our vector of trees.
  vector<Point>* centers = results->centers.get();
  vector<int>* clusterAssign = results->clusterAssignment.get();
  vector<Point>* dataPoints = data.get();
  int numCenters = centers->size();

  //cout << "There are " << numCenters << " centers\n";
  //printKMeansResults(results.get());
  string phylipName = "temp.ph";
  bool success = false;

  vector<RootedGuideTree*> subTrees;
  subTrees.resize(numCenters, 0);

  /**
   * Create a 2D vector containing all the clusters. Each row in the vector is
   * one of the k clusters.
   */
  vector<vector<int> > clusters;
  clusters.resize(numCenters);

  for(int i = 0; i < clusters.size(); i++)
    {
      //cout << "CLUSTER " << i << ": ";
      for(int j = 0; j < clusterAssign->size(); j++)
        {
	  if((*clusterAssign)[j] == i)
            {
	      //cout << j << " ";
	      clusters[i].push_back(j);
            }
        }
      //cout << "\n";
    }

  string seq = "";
  string title = "";
  string name = "";
  /**
   * Do UPGMA on each of the clusters
   */
  for(int i = 0; i < clusters.size(); i++)
    {
      cout << "Doing Cluster " << i << " Size = ";
      int sizeCluster = clusters[i].size();
      cout << sizeCluster << "\n";

      /**
       * We only perform UPGMA if there is more than 1 sequence in the group.
       */
      if(sizeCluster > 1)
        {
	  Alignment tempAlign;
	  /**
	   * Step1 create a vector<Sequence>, should only have names.
	   */
	  vector<Sequence> tempSeqVec;
	  tempSeqVec.resize(sizeCluster, Sequence(seq, name, title));
	  //cout << "The number of Seqs is " << alignObj->getNumSeqs() << "\n";
	  for(int j = 0; j < sizeCluster; j++)
            {
	      name = alignObj->getName(clusters[i][j] + 1);
	      //cout << j << ") " << name << "\n";
	      //if(name == "")
	      //{
	      //    cout << "ASDadsadsasdfsadf\n\n";
	      //}
	      tempSeqVec[j] = Sequence(seq, name, title);
            }

	  /**
	   * Step 2: create the alignment object using this.
	   */
	  tempAlign.addSequences(&tempSeqVec);
	  tempSeqVec.clear();

	  /**
	   * Step 3: Create a dist matrix from a euclidean distance of the data
               (seedMap vectors).
	  */
	  DistMatrix tempDistMat;
	  tempDistMat.ResizeRect(sizeCluster + 1);

	  // Calculate dist matrix from data, remember it is from 1 to sizeCluster
	  for (int si = 0; si < sizeCluster; si++)
            {
	      for (int sj = si + 1; sj < sizeCluster; sj++)
                {
		  tempDistMat(si + 1, sj + 1) =
		    euclideanDistance(&(*dataPoints)[clusters[i][si]],
				      &(*dataPoints)[clusters[i][sj]]);
                }
            }
	  /**
	   * Step 4: Make a call to generateTreeFromDistMatUPGMA
	   */
	  auto_ptr<RootedGuideTree> guideTree;
	  guideTree.reset(new RootedGuideTree);
	  generateTreeFromDistMatUPGMA(guideTree.get(), &tempDistMat, &tempAlign, 1,
				       sizeCluster, &phylipName, &success, false,
				       &clusters[i]);

	  subTrees[i] = guideTree.get(); // Store the subtree;
	  guideTree.release();
	  if(success)
            {
	      cout << "it was successful!\n";
            }
        }
      else
        {
	  /**
	   * There is only one sequence in the cluster. Just put a node into the position.
	   */
	  cout << "only one node\n";
	  Node* root = new Node(clusters[i][0], 0, 0);
	  root->realIndex = clusters[i][0] + 1;
	  //root
	  subTrees[i] = new RootedGuideTree(root);
        }
    }

  /**
   * At this point I have calculated all of the sub trees. I now need to do UPGMA on the
   * overall tree using the centers. I will then join everything together using this tree.
   */

  if(numCenters > 1)
    {
      Alignment tempAlign;
      name = "";
      /**
       * Step1 create a vector<Sequence>, should only have names.
       */
      vector<Sequence> tempSeqVec;
      tempSeqVec.resize(numCenters, Sequence(seq, name, title));

      for(int j = 0; j < numCenters; j++)
        {
	  ostringstream namess;
	  namess << "CLUSTER" << j;
	  name = namess.str();
	  tempSeqVec[j] = Sequence(seq, name, title);
        }

      /**
       * Step 2: create the alignment object using this.
       */
      tempAlign.addSequences(&tempSeqVec);
      tempSeqVec.clear();
      //tempAlign.printInfo();
      /**
       * Step 3: create a distance matrix from the centers.
       */
      DistMatrix centersDMat;
      centersDMat.ResizeRect(numCenters + 1);

      for (int si = 0; si < numCenters; si++)
        {
	  for (int sj = si + 1; sj < numCenters; sj++)
            {
	      centersDMat(si + 1, sj + 1) = euclideanDistance(&(*centers)[si],
							      &(*centers)[sj]);
            }
        }
      /**
       * Step 4: Generate the guide tree using UPGMA
       */
      auto_ptr<RootedGuideTree> guideTree;
      guideTree.reset(new RootedGuideTree);
      auto_ptr<AlignmentSteps> steps = generateTreeFromDistMatUPGMA(guideTree.get(),
								    &centersDMat, &tempAlign, 1,
								    numCenters, &phylipName, &success, false);
      RootedGuideTree* temp = guideTree.release();
      delete temp;
      //steps->printAlignSteps();

      /**
       * Step 5: Build up the whole guide tree from the subtrees using the new tree.
       */
      string step;
      for(int stepi = 1; stepi <= steps->getNumSteps(); stepi++)
        {
	  vector<int>* curStep = steps->getStep(stepi);
	  int indexI = -1;
	  int indexJ = -1;

	  for(int nodei = 1; nodei < curStep->size(); nodei++)
            {
	      if(indexI != -1 && indexJ != -1)
                {
		  break;
                }
	      if((*curStep)[nodei] == 1 && indexI == -1)
                {
		  indexI = nodei - 1; // This is where the left tree will be.
                }
	      if((*curStep)[nodei] == 2 && indexJ == -1)
                {
		  indexJ = nodei - 1; // This is where the right tree will be.
                }
            }
	  //cout << "Joining subcluster " << indexI << " to " << indexJ << "\n";
	  // Now indexI and indexJ should have the indices of the subtrees to
	  // be joined.
	  Node* newRoot = new Node(-1, 0, -1);
	  newRoot->left = subTrees[indexI]->getRoot();
	  newRoot->right = subTrees[indexJ]->getRoot();
	  //cout << "LEFT " << newRoot->left->realIndex << "\n";
	  //cout << "RIGHT " << newRoot->right->realIndex << "\n";
	  subTrees[indexJ] = 0;
	  subTrees[indexI] = new RootedGuideTree(newRoot);
        }
      guideTree.reset(subTrees[0]);
      printTree(guideTree.get(), 1, alignObj);
      return guideTree;
    }
  else
    {
      /**
       * Return a tree with the single node as root.
       */
      auto_ptr<RootedGuideTree> overAllGuideTree;
      overAllGuideTree.reset(subTrees[0]);
      subTrees[0] = 0;
      return overAllGuideTree;
    }
  auto_ptr<RootedGuideTree> overAllGuideTree;
  return overAllGuideTree;
}

void TreeInterface::printTree(RootedGuideTree* tree, int fSeq, Alignment* align)
{
  SeqInfo info;
  info.firstSeq = 1;
  info.lastSeq = align->getNumSeqs();
  info.numSeqs = info.lastSeq;
  RootedTreeOutput myTreeOutput(&info);

  OutputFile phylipPhyTreeFile;
  string name = "kmeans.tree";
  string path = "";

  phylipPhyTreeFile.openFile(&name, "\nEnter name for new GUIDE TREE           file  ",
                             &path, "dnd", "Guide tree");

  ofstream* ptrToFile = phylipPhyTreeFile.getPtrToFile();
  myTreeOutput.printPhylipTreeNoDists(tree, ptrToFile, align);
}

/* <-- FS */

/* FS, 2009-04-22 --> */
auto_ptr<AlignmentSteps> 
TreeInterface::generateTreeFromDistMatUPGMA(RootedGuideTree* guideTree, DistMatrix* distMat, 
					    Alignment *alignPtr, int seq1, int nSeqs,
					    string* phylipName, bool* success, bool printTree,
					    vector<int>* realNums)
{
  auto_ptr<AlignmentSteps> progSteps;
  string copyOfPhylipName = string(*phylipName);
  
  if (nSeqs >= 2)
    {
      RootedClusterTree clusterTree;
      progSteps = clusterTree.treeFromDistMatrix(guideTree, distMat, alignPtr, seq1,
						 nSeqs, copyOfPhylipName, realNums, printTree);
      
      *phylipName = copyOfPhylipName;
      cout << "Guide tree file created:   [" <<
	phylipName->c_str() << "]\n";
    }
  *success = true;
  return progSteps;
}
/* <-- FS */
