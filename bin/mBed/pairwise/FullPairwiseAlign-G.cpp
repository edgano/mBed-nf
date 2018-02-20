/* 
 * File: FullPairwiseAlign-G.cpp
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d after the last function
 *       in FullPairwiseAlign-G.cpp but inside the clustalw namespace{} 
 *
 *       clustalw namespace{
 *         firstFunction(){}
 *           ...
 *         lastFunction(){}
 *       #include "FullPairwiseAlign-G.cpp"
 *       }
 * Date: 2009-06-17 (FS)
 */


/* FS, 2009-03-30 --> */
void FullPairwiseAlign::pairwiseAlign(Alignment *alignPtr, double *dist, int iI, int iJ)
{
    int si, sj, i;
    int n, m, len1, len2;
    int maxRes;
    int _matAvgScore;
    int res;
    double _score;
    float gapOpenScale, gapExtendScale;
    
    try
    {
      if((iI < 0) || (iJ < 0))
        {
            cout << "The range for pairwise Alignment is incorrect.\n"
                 << "Need to terminate program.\n";
            exit(1);
        }
        
        _maxAlnLength = alignPtr->getMaxAlnLength();
    
        int _numSeqs = alignPtr->getNumSeqs();
        if(_numSeqs == 0)
        {
            return;
        }
    
        int num = (2 * _maxAlnLength) + 1;
        bool _DNAFlag = userParameters->getDNAFlag();
        float _pwGapOpen, _pwGapExtend;
        _pwGapOpen = userParameters->getPWGapOpen();
        _pwGapExtend = userParameters->getPWGapExtend();
        
        displ.resize(num);
        HH.resize(_maxAlnLength);
        DD.resize(_maxAlnLength);
        RR.resize(_maxAlnLength);
        SS.resize(_maxAlnLength);
        // Note these 2 lines replace the stuff above because it is all done in the SubMatrix
        PairScaleValues scaleValues;
        maxRes = subMatrix->getPairwiseMatrix(matrix, scaleValues, _matAvgScore);
        if (maxRes == 0)
        {
            cout << "Could not get the substitution matrix\n";
            return;
        }
        
        intScale = scaleValues.intScale;
        gapOpenScale = scaleValues.gapOpenScale;
        gapExtendScale = scaleValues.gapExtendScale;
    
        int _gapPos1, _gapPos2;
        _gapPos1 = userParameters->getGapPos1();
        _gapPos2 = userParameters->getGapPos2();
        const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 
    
        /*for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)*/
	si = iI;
        {
            n = alignPtr->getSeqLength(si + 1);
            len1 = 0;
            for (i = 1; i <= n; i++)
            {
                res = (*_ptrToSeqArray)[si + 1][i];
                if ((res != _gapPos1) && (res != _gapPos2))
                {
                    len1++;
                }
            }

            /*for (sj = utilityObject->MAX(si+1, jStart+1); sj < _numSeqs && sj < jEnd; sj++)*/
	    sj = iJ;
            {
                m = alignPtr->getSeqLength(sj + 1);
                if (n == 0 || m == 0)
                {
		  /*distMat->SetAt(si + 1, sj + 1, 1.0);
		    distMat->SetAt(sj + 1, si + 1, 1.0);
		    continue;*/
		  *dist = 1.0;
		  /* must free memory, FS, 2009-03-30 */
		  return; 
                }
                len2 = 0;
                for (i = 1; i <= m; i++)
                {
                    res = (*_ptrToSeqArray)[sj + 1][i];
                    if ((res != _gapPos1) && (res != _gapPos2))
                    {
                        len2++;
                    }
                }

                if (_DNAFlag)
                {
                    _gapOpen = static_cast<int>(2 * _pwGapOpen * intScale *
                                    gapOpenScale);
                    _gapExtend = static_cast<int>(_pwGapExtend * intScale * gapExtendScale);
                }
                else
                {
                    if (_matAvgScore <= 0)
                    {
                        _gapOpen = 2 * static_cast<int>((_pwGapOpen +
                               log(static_cast<double>(utilityObject->MIN(n, m)))) * intScale);
                    }
                    else
                    {
                        _gapOpen = static_cast<int>(2 * _matAvgScore * (_pwGapOpen +
                        log(static_cast<double>(utilityObject->MIN(n, m)))) * gapOpenScale);
                    }
                    _gapExtend = static_cast<int>(_pwGapExtend * intScale);
                }
                // align the sequences
            
                seq1 = si + 1;
                seq2 = sj + 1;

                _ptrToSeq1 = alignPtr->getSequence(seq1);
                _ptrToSeq2 = alignPtr->getSequence(seq2);
            
                forwardPass(_ptrToSeq1, _ptrToSeq2, n, m);
                reversePass(_ptrToSeq1, _ptrToSeq2);

                lastPrint = 0;
                printPtr = 1;

                // use Myers and Miller to align two sequences 

                maxScore = diff(sb1 - 1, sb2 - 1, se1 - sb1 + 1, se2 - sb2 + 1,
                    (int)0, (int)0);

                // calculate percentage residue identity

                mmScore = tracePath(sb1, sb2);

                if (len1 == 0 || len2 == 0)
                {
                    mmScore = 0;
                }
                else
                {
                    mmScore /= (float)utilityObject->MIN(len1, len2);
                }

                _score = ((float)100.0 - mmScore) / (float)100.0;
		*dist = _score;
                /* distMat->SetAt(si + 1, sj + 1, _score);
		   distMat->SetAt(sj + 1, si + 1, _score);*/
                
                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d",
                                        si+1, sj+1, (int)mmScore);     
                }
            }
        }

        displ.clear();
        HH.clear();
        DD.clear();
        RR.clear();
        SS.clear();
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}
/* <-- FS */
