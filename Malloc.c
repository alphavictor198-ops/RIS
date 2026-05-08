#include <stdio.h>
#include <string.h>

// Functions for firstFit, bestFit, and worstFit would be defined here
// to process the blocks and process sizes, printing the allocation results.

void firstFit(int blockSize[], int m, int processSize[], int n);
void bestFit(int blockSize[], int m, int processSize[], int n);
void worstFit(int blockSize[], int m, int processSize[], int n);

int main() {
    int bSize[] = {100, 500, 200, 300, 600};
    int pSize[] = {212, 417, 112, 426};
    int m = sizeof(bSize) / sizeof(bSize[0]);
    int n = sizeof(pSize) / sizeof(pSize[0]);

    // Create copies of block sizes to maintain original data for each algorithm
    int bSize1[m], bSize2[m], bSize3[m];
    for(int i=0; i<m; i++) bSize1[i] = bSize2[i] = bSize3[i] = bSize[i];

    firstFit(bSize1, m, pSize, n);
    bestFit(bSize2, m, pSize, n);
    worstFit(bSize3, m, pSize, n);

    return 0;
}
