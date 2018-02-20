/* 
 * File: FastPairwiseAlign-G.cpp
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d after the last function
 *       in FastPairwiseAlign-G.cpp but inside the clustalw namespace{} 
 *
 *       clustalw namespace{
 *         firstFunction(){}
 *           ...
 *         lastFunction(){}
 *       #include "FastPairwiseAlign-G.cpp"
 *       }
 * Date: 2009-06-17 (FS)
 */



/* FS, 2009-03-30 --> */
void FastPairwiseAlign::pairwiseAlign(Alignment *alignPtr, double *dist, int iI, int iJ)
{
    try
    {
        if((iI < 0) || (iJ < 0))
        {
            cout << "The range for pairwise Alignment is incorrect.\n"
                 << "Need to terminate program.\n";
            exit(1);
        }
    
        int i, j, dsr;
        double calcScore;
        bool _DNAFlag = userParameters->getDNAFlag();
        _maxAlnLength = alignPtr->getMaxAlnLength();
        int num = (2 * _maxAlnLength) + 1;
        accum.ResizeRect(5, num);
    
        displ.resize(num);
        slopes.resize(num);
        diagIndex.resize(num);

        zza.resize(_maxAlnLength + 1);
        zzb.resize(_maxAlnLength + 1);
        zzc.resize(_maxAlnLength + 1);
        zzd.resize(_maxAlnLength + 1);
    
        if (_DNAFlag)
        {
            userParameters->setDNAParams();
        }
        else
        {
            userParameters->setProtParams();
        }

        //cout << "\n\n";(gordon fri 23 April)
    
            const vector<int>* _seqIPtr = alignPtr->getSequence(iI);
            int _seqILength = alignPtr->getSeqLength(iI);
            if (_DNAFlag)
            {
                makeNPtrs(zza, zzc, _seqIPtr, _seqILength);
            }
            else
            {
                makePPtrs(zza, zzc, _seqIPtr, _seqILength);
            }
            double _score;

                const vector<int>* _seqJPtr = alignPtr->getSequence(iJ);
                int _seqJLength = alignPtr->getSeqLength(iJ);
                if (_DNAFlag)
                {
                    makeNPtrs(zzb, zzd, _seqJPtr, _seqJLength);
                }
                else
                {
                    makePPtrs(zzb, zzd, _seqJPtr, _seqJLength);
                }
                pairAlign(_seqIPtr, _seqILength, _seqJLength);
                if (!maxSoFar)
                {
                    calcScore = 0.0;
                }
                else
                {
                    calcScore = (double)accum[0][maxSoFar];
                    if (userParameters->getPercent())
                    {
                        dsr = (_seqILength < _seqJLength) ? _seqILength
                            : _seqJLength;
                        calcScore = (calcScore / (double)dsr) *100.0;
                    }
                }
                _score = (100.0 - calcScore) / 100.0;

		*dist = _score;

                if(userParameters->getDisplayInfo())
                {
                    if (calcScore > 0.1)
                    {
                        utilityObject->info("Sequences (%d:%d) Aligned. Score: %lg",
                                            iI, iJ, calcScore);     
                    }
                    else
                    {
                        utilityObject->info("Sequences (%d:%d) Not Aligned", iI, iJ);
                    }
                }

        accum.clearArray();    
        displ.clear();
        slopes.clear();
        diagIndex.clear();

        zza.clear();
        zzb.clear();
        zzc.clear();
        zzd.clear();
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FastPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }    
}
/* <-- FS */

