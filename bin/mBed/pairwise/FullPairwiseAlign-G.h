/* 
 * File: FullPairwiseAlign-G.h
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d in the last statement 
 *       of the FullPairwiseAlign class of FullPairwiseAlign.h, 
 *       that is, inside the clustalw namespace{} and the class FullPairwiseAlign{}
 *
 *       clustalw namespace{
 *       class FullPairwiseAlign{
 *          public:
 *           ...
 *          private:
 *           ...
 *       #include "FullPairwiseAlign-G.h"
 *       };
 *       }
 * Date: 2009-06-17 (FS)
 */

public:
        void pairwiseAlign(Alignment *alignPtr, double *dist, int i, int j); /* FS, 2009-03-30 */
