


auto_ptr<AlignmentSteps> 
RootedClusterTree::treeFromDistMatrix(RootedGuideTree* phyloTree,DistMatrix* distMat, Alignment *alignPtr, 
				      int seq1, int nSeqs, string& phylipName,
				      vector<int>* realNums, bool printTree)
{
    OutputFile phylipPhyTreeFile;
    auto_ptr<AlignmentSteps> progSteps;
    try
    {
        // Test to see if the inputs are valid
        if(seq1 < 1 || nSeqs < 1)
        {
            cerr << "Invalid inputs into treeFromDistMatrix \n"
                 << "seq1 = " << seq1 << " nSeqs = " << nSeqs << "\n"
                 << "Need to end program!\n";
            exit(1); 
            return progSteps;
        }

        int i;
        float dist;
        string path;
        verbose = false;
        firstSeq = seq1;
        lastSeq = firstSeq + nSeqs - 1;
    
        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = nSeqs;
        
        phylipName = "temp.dnd";
        
        //utilityObject->getPath(userParameters->getSeqName(), &path);
        
        if(nSeqs >= 2)
        {
            string name = phylipName;
            if(!phylipPhyTreeFile.openFile(&name, 
                             "\nEnter name for new GUIDE TREE           file  ", &path, "dnd",
                             "Guide tree"))
            {
                return progSteps;
            }
            phylipName = name;                    
        }
        else
        {
            return progSteps;
        }
                
        RootedTreeOutput outputTree(&info);
        
        ofstream* ptrToFile = phylipPhyTreeFile.getPtrToFile();
        
        if (nSeqs == 2)
        {
            dist = (*distMat)(firstSeq, firstSeq + 1) / 2.0;
            if(ptrToFile->is_open())
            {
                (*ptrToFile) <<  "(" << alignPtr->getName(firstSeq) << ":" 
                             << setprecision(5)
                             << dist << "," << alignPtr->getName(firstSeq + 1) << ":" 
                             << setprecision(5) << dist <<");\n";
            }
            /**
             * Create the progSteps here.
             */
            progSteps.reset(new AlignmentSteps);
            vector<int> groups;
            groups.resize(3, 0); // Set up the step
            groups[1] = 1;
            groups[2] = 2;
            progSteps->saveSet(&groups);
            bool useRealIndices = false; 
            if(realNums != 0)
            {
                if(realNums->size() == nSeqs)
                {
                    useRealIndices = true;
                }
            }
            Node* root = new Node(-1, 0, 0);
            Node* left = new Node(2, 0, 1);  
            Node* right = new Node(3, 0, 2);
            if(useRealIndices)
            {
                left->realIndex = (*realNums)[0] + 1;
                right->realIndex = (*realNums)[1] + 1;
            }
            root->left = left;
            root->right = right;
            phyloTree->setRoot(root);                     
        }
        else
        {
            UPGMAAlgorithm clusAlgorithm;
            progSteps = clusAlgorithm.generateTree(phyloTree, distMat, &info, false, realNums);
            if(printTree)
            {
                outputTree.printPhylipTree(phyloTree, ptrToFile, alignPtr, distMat);
            }
        }
        return progSteps;
    }
    catch(const exception &ex)
    {
        cout << "ERROR: Error has occured in treeFromDistMatrix. " 
             << "Need to terminate program.\n"
             << ex.what();
        exit(1);   
    }
    catch(...)
    {
        cout << "ERROR: Error has occured in treeFromDistMatrix. " 
             << "Need to terminate program.\n";      
        exit(1);  
    }
}


