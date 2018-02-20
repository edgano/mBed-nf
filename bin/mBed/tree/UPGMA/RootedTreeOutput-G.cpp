

/* FS, 2009-04-22 --> */
void RootedTreeOutput::printPhylipTreeNoDists(RootedGuideTree* tree, ofstream* ptrToFile,
					      Alignment *alignPtr)
{
  if(!ptrToFile || !ptrToFile->is_open())
    {
      return;
    }

  // If we have only 2 sequences, use the distances in the distMat
  if (lastSeq - firstSeq + 1 == 2)
    {
      (*ptrToFile) << "(" << alignPtr->getName(firstSeq) << ","
		   << alignPtr->getName(firstSeq + 1);
      return ;
    }

  //(*ptrToFile) << "(\n";
  phylipTraverseNoHeights(ptrToFile, alignPtr, tree->getRoot());
  (*ptrToFile) << ")\n";
}
/* <-- FS */


/* FS, 2009-04-22 --> */
void RootedTreeOutput::phylipTraverseNoHeights(ofstream* ptrToFile, Alignment *alignPtr, Node* t)
{
  if(!ptrToFile)
    {
      return;
    }
  if(t != 0)
    {
      if(isLeafNode(t))
        {
	  if(alignPtr)
            {
	      (*ptrToFile) << alignPtr->getName(t->realIndex);
            }
	  else
            {
	      (*ptrToFile) << t->realIndex;
            }
        }
      else // Internal node
        {
	  (*ptrToFile) << "(\n";
	  phylipTraverseNoHeights(ptrToFile, alignPtr, t->left);
	  (*ptrToFile) << ",\n";
	  phylipTraverseNoHeights(ptrToFile, alignPtr, t->right);
	  (*ptrToFile) << ")";
        }
    }
}

bool RootedTreeOutput::isLeafNode(Node* t)
{
  if(t->left == 0 && t->right == 0)
    {
      return true;
    }
  else
    {
      return false;
    }
}

/* <-- FS */
