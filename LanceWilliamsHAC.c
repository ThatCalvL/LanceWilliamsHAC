///////////////////////////////////////////////////////////////////////

// COMP2521 Assignment 2 Part 3
// Written by Calvin Li z5242094
// 15 Nov 2020

/////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "LanceWilliamsHAC.h"


#define INFINITY INT_MAX
#define COMPLETE 2
#define SINGLE 1

/////////////////////////////////// Helper functions //////////////////////////////////////
// find the smaller number between a and b
static double min(double a,double b);

// find the larger number between a and b
static double max(double a, double b);

// create a newDNode to join cluster A and cluster B
static Dendrogram newDNode (int a, int b, Dendrogram *Dgram);

//given methods and clusters (with incremental interger and distance array), update according to method
static void methodUpdate(int method, int placeA, int placeB, int incre, double** distance);

// merge valid DNode 
static void mergeCluster (Graph g, Dendrogram *Dgram, int placeA, int placeB, int method, double** distance);

// find shortest distance and merge DNode
static int findShortestAndMerge(Graph g, double **distance, Dendrogram *Dgram, int method);



Dendrogram LanceWilliamsHAC(Graph g, int method){

    // Create a matrix for recording dist[a,b]
    double** distance = malloc(GraphNumVertices(g)*sizeof(double*));
    int j;
    for (j = 0; j < GraphNumVertices(g); j++){
        distance[j] = malloc(GraphNumVertices(g)*sizeof(double));
    }

    // Create a array of Dendrogram Node
    Dendrogram *Dgram = malloc(GraphNumVertices(g)*sizeof(Dendrogram));
    int i;
    for (i = 0; i < GraphNumVertices(g); i++){
        Dgram[i] = malloc(sizeof(DNode));
        Dgram[i]->vertex = i;
        Dgram[i]->left = NULL;
        Dgram[i]->right = NULL;
    }

    // create a matrix to record direct distance between 2 vertex.
    // set all value in direct matrix to 0;
    int** direct = malloc(GraphNumVertices(g)*sizeof(int*));
    int ss;
    for (ss = 0; ss < GraphNumVertices(g); ss++){
        direct[ss] = malloc(GraphNumVertices(g)*sizeof(int));
        int kk;
        for (kk = 0; kk < GraphNumVertices(g); kk++){
            direct[ss][kk] = 0;
        }   
    }
    
    // find and fill in direct matrix with edges between two vertices.
    int r;
    for (r = 0; r < GraphNumVertices(g); r++){
        AdjList record = GraphOutIncident(g,r);
        while(record!= NULL){
            direct[r][record->v] = record->weight;
            record = record->next;
        }
    }
   
    //record the index of DNode to return
    int recordLast = 0;

        int z;
        for (z = 0; z < GraphNumVertices(g); z++){
        // create ladder matrix to get 1/max(weight) of direct edges
            int m;
            for (m = z+1; m < GraphNumVertices(g); m++){
                if(direct[z][m] == 0 && direct[m][z] == 0){
                // if no direct edges
                    distance[z][m] = INFINITY;
                } else{
                    double dA = (double)direct[m][z];
                    double dB = (double)direct[z][m];
                    distance[z][m] = 1/max(dA,dB);
                }
            }
        }
        int p = findShortestAndMerge(g,distance,Dgram,method);
        
        // start to find cluster and merge DNnode
        while(p != -1){
            // constantly find clusters and merge until no cluster to be merged
            recordLast = p;
            p = findShortestAndMerge(g,distance,Dgram,method);
        }

    return Dgram[recordLast];
    // return the last DNode 
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
    free(d);
}

/////////////////////////////////// Helper implementation //////////////////////////////////////
// for single link method method
// return the smaller one
static double min(double a,double b){
    if(a == INFINITY && b == INFINITY) return INFINITY; 
    if(a == INFINITY) return b;
    if(b == INFINITY) return a;
    if(a > b) return b;
    return a;
}
// for completed link method
// return the bigger one
static double max(double a, double b){
    if(a == INFINITY && b == INFINITY) return INFINITY;
    if(a == INFINITY) return b;
    if(b == INFINITY) return a;
    if(a > b) return a;
    return b;
}

static Dendrogram newDNode (int a, int b, Dendrogram *Dgram) {
    Dendrogram newGram = malloc(sizeof(DNode));
    newGram->vertex = -1;
    newGram->left = Dgram[a];
    newGram->right = Dgram[b];
    
    return newGram;
}

// find shortest distance and merge DNode
static int findShortestAndMerge(Graph g, double **distance, Dendrogram *Dgram, int method){
    // initialize singleLink 
    double shortest = INFINITY;
    int placeA = -1;
    int placeB = -1;

    // find the shortest path between clusters
    for (int k = 0; k < GraphNumVertices(g); k++){
        for (int m = k+1; m < GraphNumVertices(g); m++){
            if(distance[k][m] == 0){
                continue;
            } else if(distance[k][m] <= shortest){
                shortest = distance[k][m];
                placeA = k;
                placeB = m;
            }
        }
    }

    // MERGE PRECESS
    if(placeA == -1 && placeB == -1){
    // all the clusters has been merged
        return -1;
    } else{
        mergeCluster (g, Dgram, placeA, placeB, method, distance);
    }
    return placeA;
}

// merge valid DNode 
static void mergeCluster (Graph g, Dendrogram *Dgram, int placeA, int placeB, int method, double** distance) {
    // create a newDNode to join cluster A and cluster B
    Dendrogram newGram = newDNode(placeA, placeB, Dgram);
    //update the array
    //put the new Dgram in the smaller index vertex
    Dgram[placeA] = newGram;
    Dgram[placeB] = NULL;

    //update the matrix
    int q = 0;
    while(q < GraphNumVertices(g)){
        // since the array is a ladder so there are several conditions for q
        if(q == placeA){
            // renew distance[A][B] to be set as found
            distance[q][placeB] = 0;
        }
        methodUpdate(method, placeA, placeB, q, distance);
        q++;
    }
}

//given methods and clusters (with incremental interger and distance array), update according to method
static void methodUpdate(int method, int placeA, int placeB, int incre, double** distance) {
    if(incre < placeA){
        if(method == COMPLETE){
            distance[incre][placeA] = max(distance[incre][placeA],distance[incre][placeB]);
        }
        if(method == SINGLE){
            distance[incre][placeA] = min(distance[incre][placeA],distance[incre][placeB]);
        }
        distance[incre][placeB] = 0;
        // find second cluster and merge it to the smaller vertex A
    }
    if(incre > placeA && incre < placeB){
        if(method == COMPLETE){
            distance[placeA][incre] = max(distance[placeA][incre],distance[incre][placeB]);
        }
        if(method == SINGLE){
            distance[placeA][incre] = min(distance[placeA][incre],distance[incre][placeB]);
        }
        distance[incre][placeB] = 0;
    }
    if(incre == placeB){ //do nothing
    }
    if(incre > placeB){
        if(method == COMPLETE){
            distance[placeA][incre] = max(distance[placeA][incre],distance[placeB][incre]);
        }
        if(method == SINGLE){
            distance[placeA][incre] = min(distance[placeA][incre],distance[placeB][incre]);
        }
        distance[placeB][incre] = 0;
    }
}
