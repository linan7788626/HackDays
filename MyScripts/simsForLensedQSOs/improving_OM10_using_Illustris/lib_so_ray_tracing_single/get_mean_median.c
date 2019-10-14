//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

double GetMedian(double *daArray, int iSize) {
    // Allocate an array of the same size and sort it.
	int i,j;
    double * dpSorted = (double *)malloc(iSize*sizeof(double));
    for (i = 0; i < iSize; ++i) {
        dpSorted[i] = daArray[i];
    }
	double dTemp = 0.0;
    for (i = iSize - 1; i > 0; --i) {
        for (j = 0; j < i; ++j) {
            if (dpSorted[j] > dpSorted[j+1]) {
                dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    double dMedian = 0.0;
    if ((iSize % 2) == 0) {
        dMedian = (dpSorted[iSize/2] + dpSorted[(iSize/2) - 1])/2.0;
    } else {
        dMedian = dpSorted[iSize/2];
    }
    free(dpSorted);
    return dMedian;
}

double GetMode(double *daArray, int iSize) {
    // Allocate an int array of the same size to hold the
    // repetition count
	int i,j;
    int* ipRepetition = (int *)malloc(iSize*sizeof(int));
    for (i = 0; i < iSize; ++i) {
        ipRepetition[i] = 0;
        j = 0;
        //bool bFound = false;
        while ((j < i) && (daArray[i] != daArray[j])) {
            if (daArray[i] != daArray[j]) {
                ++j;
            }
        }
        ++(ipRepetition[j]);
    }
    int iMaxRepeat = 0;
    for (i = 1; i < iSize; ++i) {
        if (ipRepetition[i] > ipRepetition[iMaxRepeat]) {
            iMaxRepeat = i;
        }
    }
    free(ipRepetition);
    return daArray[iMaxRepeat];
}

double GetMean(double *daArray, int iSize) {
    double dSum = daArray[0];
	int i;
    for (i = 1; i < iSize; ++i) {
        dSum += daArray[i];
    }
    return dSum/iSize;
}
