/* 
 * File: FastPairwiseAlign-G.h
 *
 * Note: this file contains changes necessary to the original code
 * Note: this file should be #include-d in the last statement 
 *       of the FastPairwiseAlign class of FastPairwiseAlign.h, 
 *       that is, inside the clustalw namespace{} and the class FastPairwiseAlign{}
 *
 *       clustalw namespace{
 *       class TreeInterface{
 *          public:
 *           ...
 *          private:
 *           ...
 *       #include "FastPairwiseAlign-G.h"
 *       };
 *       }
 * Date: 2009-06-17 (FS)
 */


public:
    void pairwiseAlign(Alignment *alignPtr, double *dist, int i, int j); /* FS, 2009-03-30 */
