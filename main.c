//
//  main.c
//  CubicLat
//
//  Created by Simone Longo on 10/1/18.
//  Copyright © 2018 Simone Longo. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Declare global variables */

gsl_rng *rng; // Random number generator
FILE *xyzFP; // xyz file with specified number of frames
FILE *clusterFP;
FILE *outputFP;

int ****pos;
int **parts;

double CovE, VdwE;
double VdwScale = 0.25;
double defaultCovE = -4.5;
double defaultVdwE = -1.0;
int defaultNumSweeps = 1000000;

int size, npartLink, npartMod, npart;

// The desired number of output frames spread evenly throughout the simulation
int NUM_FRAMES = 1000;

// Use these to refer to properties specific to certain particles
// Linker has 3 rotation states: (0) aligned with x-axis, (1) in y, (2) in z
// Central molecules have only 1 rotation state, i.e. aligned with all axes.
enum partType {
    EMPTY, LINKER, MODULATOR
};
enum partIndex {
    XCOORD, YCOORD, ZCOORD, TYPE, ROTATION, VISITED
};
enum posIndex {
    POSTYPE, POSROT, POSVISITED
};

// Global variables used in cluster analysis
int ***clusterList;
int **tempCluster;

// Prototype for output file
void printXYZ(void);

void printArgumentError(int, int);

// Next few lines taken from https://stackoverflow.com/questions/4217037/catch-ctrl-c-in-c
// This method catches termination signals and closes the file so data is not lost
static volatile int keepRunning = 1;

void intHandler(int dummy) {
    keepRunning = 0;
    printf("Program exited due to control %d\n", dummy);
    fflush(stdout);
    fflush(clusterFP);
    fflush(xyzFP);
    fflush(outputFP);
    fclose(outputFP);
    fclose(clusterFP);
    fclose(xyzFP); // Close file
    exit(0);
}

/* Initiates pos to zero*/
void resetPos() {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                for (int l = 0; l < 3; l++) {
                    pos[i][j][k][l] = 0;
                }
            }
        }
    }
}

/* Uses specified number of particles and randomly assigns positions */
void molInit() {
    for (int i = 0; i < npartMod; i++) {
        parts[i][TYPE] = MODULATOR;
        parts[i][ROTATION] = (int) gsl_rng_uniform_int(rng, 3);
        parts[i][VISITED] = 0;

        int x = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        int y = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        int z = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        // Make sure site isn't occupied
        while (pos[x][y][z][POSTYPE] != EMPTY) {
            x = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
            y = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
            z = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        }
        // Update pos and parts with positions
        parts[i][XCOORD] = x;
        parts[i][YCOORD] = y;
        parts[i][ZCOORD] = z;
        pos[x][y][z][POSTYPE] = parts[i][TYPE];
        pos[x][y][z][POSROT] = parts[i][ROTATION];
    }
    // Now do same for LINKER types
    for (int i = npartMod; i < npartMod + npartLink; i++) {
        parts[i][TYPE] = LINKER;
        // Since LINKER is 2D, orientation on axis is important
        parts[i][ROTATION] = (int) gsl_rng_uniform_int(rng, 6);
        parts[i][VISITED] = 0;

        int x = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        int y = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        int z = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        // Make sure site isn't occupied
        while (pos[x][y][z][POSTYPE] != EMPTY) {
            x = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
            y = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
            z = (int) gsl_rng_uniform_int(rng, (unsigned long)size);
        }
        // Update pos and parts with positions
        parts[i][XCOORD] = x;
        parts[i][YCOORD] = y;
        parts[i][ZCOORD] = z;
        pos[x][y][z][POSTYPE] = parts[i][TYPE];
        pos[x][y][z][POSROT] = parts[i][ROTATION];
    }
}


/* Returns x, y, and z coordinates of neighboring particle in a specified direction */
void getNeighborIndex(int position[3], int dir) {
    switch (dir) {
        case 0: // Move left along x-axis
            position[XCOORD] = (position[XCOORD] - 1 + size) % size;
            break;
        case 1: // Move up along y-axis
            position[YCOORD] = (position[YCOORD] + 1 + size) % size;
            break;
        case 2: // Move right along x-axis
            position[XCOORD] = (position[XCOORD] + 1 + size) % size;
            break;
        case 3: // Move down along y-axis
            position[YCOORD] = (position[YCOORD] - 1 + size) % size;
            break;
        case 4: // Move forward (out of plane) along z-axis
            position[ZCOORD] = (position[ZCOORD] + 1 + size) % size;
            break;
        case 5: // Move back along z-axis
            position[ZCOORD] = (position[ZCOORD] - 1 + size) % size;
            break;
        default:
            printf("ERROR: particle can only move in one of 6 directions. %d is an invalid direction.\n", dir);
            break;
    }
}

/* Returns 1 (true) if the 2 input particles have the appropriate orientations for bonding. Returns 0 (false) otherwise */
int bondCompatible(enum partType type1, int rot1, enum partType type2, int rot2, int dir) {
    if (type1 == EMPTY || type2 == EMPTY) {
        return 0;
    }
    if (type1 == type2) {
        return 0;
    }
    if (type1 == MODULATOR) {
        switch (dir) {
            case 0:
                if ((rot1 == 0 && rot2 == 0) || (rot1 == 1 && rot2 == 1)) return 1;
                else return 0;
            case 1:
                if ((rot1 == 0 && rot2 == 2) || (rot1 == 2 && rot2 == 3)) return 1;
                else return 0;
            case 2:
                if ((rot1 == 0 && rot2 == 0) || (rot1 == 1 && rot2 == 1)) return 1;
                else return 0;
            case 3:
                if ((rot1 == 0 && rot2 == 2) || (rot1 == 2 && rot2 == 3)) return 1;
                else return 0;
            case 4:
                if ((rot1 == 1 && rot2 == 4) || (rot1 == 2 && rot2 == 5)) return 1;
                else return 0;
            case 5:
                if ((rot1 == 1 && rot2 == 4) || (rot1 == 2 && rot2 == 5)) return 1;
                else return 0;
            default:
                printf("ERROR in bondcompatible where type1 is MODULATOR.\n");
                return 0;
        }
    }
    else { // type1 is a LINKER and type2 is MODULATOR
        switch (dir) {
            case 0:
                if ((rot1 == 0 && rot2 == 0) || (rot1 == 1 && rot2 == 1)) return 1;
                else return 0;
            case 1:
                if ((rot1 == 2 && rot2 == 0) || (rot1 == 3 && rot2 == 2)) return 1;
                else return 0;
            case 2:
                if ((rot1 == 0 && rot2 == 0) || (rot1 == 1 && rot2 == 1)) return 1;
                else return 0;
            case 3:
                if ((rot1 == 2 && rot2 == 0) || (rot1 == 3 && rot2 == 2)) return 1;
                else return 0;
            case 4:
                if ((rot1 == 4 && rot2 == 1) || (rot1 == 5 && rot2 == 2)) return 1;
                else return 0;
            case 5:
                if ((rot1 == 4 && rot2 == 1) || (rot1 == 5 && rot2 == 2)) return 1;
                else return 0;
            default:
                printf("ERROR in bondcompatible where type1 is LINKER.\n");
                return 0;
        }
    }
}

/* Returns energy of an indvidual particle based on its neighbors */
double getEnergy(int x, int y, int z, enum partType type, int rot) {
    double energy = 0.0;
    for (int i = 0; i < 6; i++) {
        int newPos[3] = {x, y, z};
        getNeighborIndex(newPos, i);

        if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] == EMPTY) {
            continue; // No particle, no party.
        }

        // If there is a particle in the position we are considering...
        switch (type) {
            case MODULATOR:
                if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] == MODULATOR) {
                    energy += VdwE;
                }
                if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] == LINKER) {
                    // If the bond is compatible
                    if (bondCompatible(type, rot, pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE],
                                       pos[newPos[0]][newPos[1]][newPos[2]][POSROT], i)) {
                        energy += CovE;
                    } else {
                        energy += VdwE * VdwScale;
                    }
                }
                break;
            case LINKER:
                if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] == LINKER) {
                    energy += VdwE;
                }
                if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] == MODULATOR) {
                    // If the bond is compatible
                    if (bondCompatible(type, rot, pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE],
                                       pos[newPos[0]][newPos[1]][newPos[2]][POSROT], i)) {
                        energy += CovE;
                    } else {
                        energy += VdwE * VdwScale;
                    }
                }
                break;
            default:
                printf("ERROR: Incorrect particle type being requested.\n");
                return 0.0;
        }
    }
    return energy;
}

/* Makes an actual move, returns the energy of that move, then moves particle back */
double getHypEnergy(int x, int y, int z, enum partType type, int rot, int dir) {
    // To pass as an argument for new coords
    int newPos[3] = {x, y, z};

    getNeighborIndex(newPos, dir);

    // Check if new position is occupied
    if (pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] != EMPTY) {
        return 0.0;
    }

    // Clear old position
    pos[x][y][z][POSTYPE] = EMPTY;
    pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] = type;
    pos[newPos[0]][newPos[1]][newPos[2]][POSROT] = rot;

    double energy = getEnergy(newPos[0], newPos[1], newPos[2], type, rot);

    // Return everything to original configuration
    pos[x][y][z][POSTYPE] = type;
    pos[x][y][z][POSROT] = rot;
    pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] = EMPTY;
    pos[newPos[0]][newPos[1]][newPos[2]][POSROT] = 0;

    return energy;
}

/* Propose a move to be accepted or rejected based on energy differences and/or random numbers */
double proposeMove(int rp) {
    int newPos[3];
    for (int i = 0; i < 3; i++) {
        // Copy x, y, z coords
        newPos[i] = parts[rp][i];
    }
    int randDirection = (int) gsl_rng_uniform_int(rng, 6); // Can move in one of 6 directions
    getNeighborIndex(newPos, randDirection);

    // If new position is occupied, ignore
    if (pos[newPos[XCOORD]][newPos[YCOORD]][newPos[ZCOORD]][POSTYPE] != EMPTY) {
        return 0.0;
    }

    double oldE = getEnergy(parts[rp][XCOORD], parts[rp][YCOORD], parts[rp][ZCOORD], parts[rp][TYPE],
                            parts[rp][ROTATION]);
    double newE = getHypEnergy(parts[rp][XCOORD], parts[rp][YCOORD], parts[rp][ZCOORD], parts[rp][TYPE],
                               parts[rp][ROTATION], randDirection);
    double dE = newE - oldE;

    if (gsl_rng_uniform(rng) <= exp(-dE)) { // Accept move
        // Clear old position
        pos[parts[rp][XCOORD]][parts[rp][YCOORD]][parts[rp][ZCOORD]][POSTYPE] = EMPTY;
        pos[parts[rp][XCOORD]][parts[rp][YCOORD]][parts[rp][ZCOORD]][POSROT] = 0;

        // Set new position
        pos[newPos[0]][newPos[1]][newPos[2]][POSTYPE] = parts[rp][TYPE];
        pos[newPos[0]][newPos[1]][newPos[2]][POSROT] = parts[rp][ROTATION];

        // Update parts
        for (int i = 0; i < 3; i++) {
            parts[rp][i] = newPos[i];
        }
        return dE;
    } else { // Reject move
        return 0.0;
    }
}

/* Propose a rotation to be accepted or rejected based on energy differences and/or random numbers */
double proposeRotation(int rp) {
    // Assign new rotation based on TYPE
    int newRot = parts[rp][TYPE] == MODULATOR ? (int) gsl_rng_uniform_int(rng, 3) : (int) gsl_rng_uniform_int(rng, 6);

    // If rotation doesn't change, neither will energy
    if (parts[rp][ROTATION] == newRot) {
        return 0.0;
    }

    double oldE = getEnergy(parts[rp][XCOORD], parts[rp][YCOORD], parts[rp][ZCOORD], parts[rp][TYPE],
                            parts[rp][ROTATION]);
    double newE = getEnergy(parts[rp][XCOORD], parts[rp][YCOORD], parts[rp][ZCOORD], parts[rp][TYPE], newRot);
    double dE = newE - oldE;

    if (gsl_rng_uniform(rng) < exp(-dE)) {
        // Accept move and update everything
        pos[parts[rp][XCOORD]][parts[rp][YCOORD]][parts[rp][ZCOORD]][POSROT] = newRot;
        parts[rp][ROTATION] = newRot;
        return dE;
    } else { // Reject move
        return 0.0;
    }
}

/* Returns total energy of system */
double getTotalEnergy() {
    double energy = 0.0;
    for (int i = 0; i < npart; i++) {
        double e = getEnergy(parts[i][XCOORD], parts[i][YCOORD], parts[i][ZCOORD], parts[i][TYPE], parts[i][ROTATION]);
        energy += e;
    }
    // Divide by 2 to account for double counting energies
    return energy / 2.0;
}

/* Returns number of covalently bound molecules */
int getNumBonds(int x, int y, int z, enum partType type, int rot) {
    int numBonds = 0;
    for (int dir = 0; dir < 6; dir++) {
        int coords[3] = {x, y, z};

        // Get new neighbor coordinates for each direction
        getNeighborIndex(coords, dir);

        // If empty, no interaction
        if (pos[coords[0]][coords[1]][coords[2]][POSTYPE] == EMPTY) {
            continue;
        }

        // If same type, no covalent interaction
        if (pos[coords[0]][coords[1]][coords[2]][POSTYPE] == type) {
            continue;
        }

        if (bondCompatible(type, rot, pos[coords[0]][coords[1]][coords[2]][POSTYPE],
                           pos[coords[0]][coords[1]][coords[2]][POSROT], dir)) {
            numBonds++;
        }
    }
    return numBonds;
}

/* Returns total number of covalent bonds in system */
int getTotalCovalent() {
    int numBonds = 0;
    for (int i = 0; i < npart; i++) {
        numBonds += getNumBonds(parts[i][XCOORD], parts[i][YCOORD], parts[i][ZCOORD], parts[i][TYPE],
                                parts[i][ROTATION]);
    }
    // Total number of bond should be even since they were all double counted
    if (numBonds % 2) {
        printf("ERROR: unequal number of covalent bonds.\n");
    }
    return numBonds / 2;
}

/*
 * -------------------------------------------------------------------------------------
 *                                    CLUSTER ANALYSIS
 * -------------------------------------------------------------------------------------
 */

/* Initialize cluster data to zero */
void initializeClusterList() {
    for (int i = 0; i < npart; i++) {
        for (int j = 0; j < npart; j++) {
            for (int k = 0; k < 6; k++) {
                clusterList[i][j][k] = 0;
            }
        }
    }
}

/* Resets all particles in temp cluster to zero */
void resetTempCluster() {
    for (int i = 0; i < npart; i++) {
        for (int j = 0; j < 6; j++) {
            tempCluster[i][j] = 0;
        }
    }
}

/* Marks a particle as visited in both parts and pos arrays */
void markVisited(int part) {
    if (part < 0) return;
    parts[part][VISITED] = 1;
    pos[parts[part][XCOORD]][parts[part][YCOORD]][parts[part][ZCOORD]][POSVISITED] = 1;
}

/* Returns the size of the input 2D array with 2nd dim containing particle representation like in parts */
int getClusterSize(int **array) {
    int i = 0;
    for (; i < npart; i++) {
        if (array[i][TYPE] == EMPTY) {
            return i;
        }
    }
    return -1;
}

/* Checks tempCluster to see if specified particle is already contained in the set.
   Returns 1 if the particle is already in the set, 0 otherwise */
int containsParticle(int part) {
    for (int i = 0; i < getClusterSize(tempCluster); i++) {
        if (parts[part][XCOORD] == tempCluster[i][XCOORD] &&
            parts[part][YCOORD] == tempCluster[i][YCOORD] &&
            parts[part][ZCOORD] == tempCluster[i][ZCOORD] &&
            parts[part][TYPE] == tempCluster[i][TYPE] &&
            parts[part][ROTATION] == tempCluster[i][ROTATION]) {
            if (parts[part][VISITED] != tempCluster[i][VISITED]) {
                printf("ERROR: Discrepant visitation!\n");
                return 0;
            }
            return 1;
        }
    }
    return 0;
}

/* Adds particle to tempCluster list AND marks as visited in parts and pos returns 1 if successful, 0 otherwise */
int addTemp(int part) {
    if (part == -1) {
        printf("\nERROR: Particle index out of bounds!\n");
        return 0;
    }
    if (containsParticle(part) && parts[part][VISITED] == 1) {
        return 0;
    }
    // This is the index where the next particle will be stored
    int tmpInd = getClusterSize(tempCluster);
    if (tmpInd == -1) {
        printf("\nERROR: Cluster too large! Line: %d\n", __LINE__);
        return 0;
    }
    markVisited(part);
    for (int i = 0; i < 6; i++) {
        tempCluster[tmpInd][i] = parts[part][i];
    }
    return 1;
}

/* Returns an array of all unvisited nearest neighbors of an individual particle */
void findNearestUnvisited(int part, int nn[][6]) {
    if (part < 0 || part > npart) {
        printf("Illegal part number!\n");
        exit(0);
    }

    // Set everything in nearest neighbors to 0
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            nn[i][j] = 0;
        }
    }

    // Check every direction for neighbors
    for (int dir = 0; dir < 6; dir++) {
        int coords[3] = {parts[part][XCOORD], parts[part][YCOORD], parts[part][ZCOORD]};
        getNeighborIndex(coords, dir);
        if (bondCompatible(parts[part][TYPE], parts[part][ROTATION], pos[coords[0]][coords[1]][coords[2]][POSTYPE],
                           pos[coords[0]][coords[1]][coords[2]][POSROT], dir)) {
            // If bond is compatible and particle has not been visited yet
            if (pos[coords[0]][coords[1]][coords[2]][POSVISITED] == 0) {
                nn[dir][XCOORD] = coords[0];
                nn[dir][YCOORD] = coords[1];
                nn[dir][ZCOORD] = coords[2];
                nn[dir][TYPE] = pos[coords[0]][coords[1]][coords[2]][POSTYPE];
                nn[dir][ROTATION] = pos[coords[0]][coords[1]][coords[2]][POSROT];
                nn[dir][VISITED] = pos[coords[0]][coords[1]][coords[2]][POSVISITED];
            }
        }
    }
}

/* Returns number of nearest unvisited neighbors */
int numUnvisitedNN(int nn[][6]) {
    int count = 0;
    for (int i = 0; i < 6; i++) {
        if (nn[i][TYPE] != EMPTY) {
            count++;
        }
    }
    return count;
}

/* Linear search through parts to return index of particle with specified coordinates. Returns -1 if not found. */
int findPartNum(int x, int y, int z) {
    for (int i = 0; i < npart; i++) {
        if (parts[i][XCOORD] == x && parts[i][YCOORD] == y && parts[i][ZCOORD] == z) {
            return i;
        }
    }
    return -1;
}

/* Returns next open index for adding a new cluster. Returns -1 if it doesn't exist. */
int nextCluster() {
    for (int i = 0; i < npart; i++) {
        if (clusterList[i][0][TYPE] == EMPTY) {
            return i;
        }
    }

    printf("ERROR: no available index for cluster.\n");
    return -1;
}

/* Recursive method called by findClusters() to find clusters local to each particle */
void findCluster(int part) {
    int nn[6][6]; // Stores nearest neighbors (if they exist) in all 6 possible directions
    findNearestUnvisited(part, nn); // Populate nearest neighbors list
    int numUnvisited = numUnvisitedNN(nn);
    // Base case
    if (numUnvisited == 0 && parts[part][VISITED] == 0) {
        addTemp(part);
        return;
    }
    for (int i = 0; i < 6; i++) {
        // If particle exists and is unvisited
        if (nn[i][TYPE] != 0 && nn[i][VISITED] == 0) {
            int nextPartInd = findPartNum(nn[i][XCOORD], nn[i][YCOORD], nn[i][ZCOORD]);
            if (nextPartInd < 0 || nextPartInd >= npart) {
                printf("ERROR: Illegal particle number!\n");
            }
            nn[i][VISITED] = 1;
            addTemp(part);
            // Recursive case
            findCluster(nextPartInd);
        }
    }
}

/* Gives coordinates of a cluster's center of mass */
void getClusterCenter(int clusterNum, double *coords) {
    int sizeCluster = getClusterSize(clusterList[clusterNum]);
    int x = 0, y = 0, z = 0;
    for (int i = 0; i < sizeCluster; i++) {
        x += clusterList[clusterNum][i][XCOORD];
        y += clusterList[clusterNum][i][YCOORD];
        z += clusterList[clusterNum][i][ZCOORD];
    }
    coords[XCOORD] = (double) x / (double) sizeCluster;
    coords[YCOORD] = (double) y / (double) sizeCluster;
    coords[ZCOORD] = (double) z / (double) sizeCluster;
}

/* Find all clustered particles with respect to all particles */
void findClusters() {
    int size1clusters = 0;
    for (int i = 0; i < npart; i++) {
        resetTempCluster();
        findCluster(i);
        int clusterSize = getClusterSize(tempCluster);
        if (clusterSize == 0) {
            continue;
        }

        // Next available index for adding a new cluster
        int addInd = nextCluster();
        for (int j = 0; j < clusterSize; j++) {
            for (int k = 0; k < 6; k++) {
                // Copy temp cluster to comprehensive cluster list
                clusterList[addInd][j][k] = tempCluster[j][k];
            }
        }
        if (clusterSize == 1) {
            size1clusters++;
            continue;
        }
        double centerCoords[3] = {-1, -1, -1};
        getClusterCenter(addInd, centerCoords);
        printf("Cluster size %d: %d\t Center: %f %f %f\n", i, clusterSize, centerCoords[0], centerCoords[1],
               centerCoords[2]);
    }

    printf("Number of unbonded particles: %d\n", size1clusters);
}

/* Print cluster information to clusterFP */
void analyzeClusters() {
    // Number of clusters
    int numClusters = nextCluster();
    int avgSize = 0;
    int maxSize = 0;
    int i = 0;
    if (numClusters == 0)
        return;

    int size1clusters = 0;
    fprintf(clusterFP, "ClustNum\tClustSize\tClustCenter\n");
    for (i = 0; i < numClusters; i++) {
        int size = getClusterSize(clusterList[i]);
        if (size <= 0) {
            continue;
        }
        if (size == 1) {
            size1clusters++;
            continue;
        }
        double centerCoords[3] = {-1, -1, -1};
        getClusterCenter(i, centerCoords);
        fprintf(clusterFP, "%d\t%d\t%f, %f, %f\n", (i + 1), size, centerCoords[0], centerCoords[1], centerCoords[2]);
        avgSize += size;
        maxSize = size > maxSize ? size : maxSize;
    }
    fprintf(clusterFP, "\nTotal number of clusters: %d\n", i);
    fprintf(clusterFP, "Average Cluster Size: %f\n", ((double) avgSize / (double) i));
    fprintf(clusterFP, "Biggest cluster: %d\n", maxSize);
    // Note: avgSize as used below is intended to show total number of clustered particles, not actual average.
    fprintf(clusterFP, "Total particles in a cluster: %d out of %d system particles.\n", avgSize, npart);
}



/*
 * -------------------------------------------------------------------------------------
 *                                  END CLUSTER ANALYSIS
 *                                Thank you for visiting...
 * -------------------------------------------------------------------------------------
 */

/*
 *
 *  MAIN METHOD
 *  ------------------------------------------
 *  Takes 6 arguments:
 *  (1) Size
 *  (2) Number of linker molecules
 *  (3) Number of central molecules
 *  (4) Covalent bond energy
 *  (5) Non specific interaction energy
 *  (6) Desired number of Monte Carlo sweeps
 *  ------------------------------------------
 */
int main(int argc, const char *argv[]) {
    /* Parsing arguments */
    int numSweeps = defaultNumSweeps;
    // If all parameters have been passed
    if (argc == 7) {
        size = atoi(argv[1]);
        npartLink = atoi(argv[2]);
        npartMod = atoi(argv[3]);
        CovE = atof(argv[4]);
        VdwE = atof(argv[5]);
        numSweeps = atoi(argv[6]);
        npart = npartLink + npartMod;
    }
    // If only 4 are passed
    if (argc == 5) {
        size = atoi(argv[1]);
        npartLink = atoi(argv[2]);
        npartMod = atoi(argv[3]);
        CovE = defaultCovE;
        VdwE = defaultVdwE;
        numSweeps = atoi(argv[4]);
        npart = npartLink + npartMod;
    }
    // If only 3 are passed
    if (argc == 4) {
        size = atoi(argv[1]);
        npartLink = atoi(argv[2]);
        npartMod = atoi(argv[3]);
        CovE = defaultCovE;
        VdwE = defaultVdwE;
        numSweeps = defaultNumSweeps;
        npart = npartLink + npartMod;
    }
    if (numSweeps < 1000) {
        printArgumentError(numSweeps, 1);
        exit(0);
    } else if (argc != 7 && argc != 5 && argc != 4) {
        printArgumentError(argc, 0);
        exit(0);
    }

    xyzFP = fopen("pos.lammpstrj", "w");
    outputFP = fopen("output.txt", "w");

    /* RNG environment setup */
    int seed = 0;
    gsl_rng_env_setup();
    /* TO USE:
     * (int) gsl_rng_uniform_int(rng, <non-inclusive max int value>)
     * (double) gsl_rng_uniform(rng)
     *
     * To compile: make sure to add -lgsl flag
     */
    rng = gsl_rng_alloc(gsl_rng_taus2);
    if (seed == 0) {
        seed = (int) time(NULL);
    }
    printf("seed = %i\n", seed);
    gsl_rng_set(rng, (unsigned long)seed);

    // Memory allocation for data structures being used
    pos = (int ****) malloc(size * sizeof(int ***));
    for (int i = 0; i < size; i++) {
        pos[i] = (int ***) malloc(size * sizeof(int **));
        for (int j = 0; j < size; j++) {
            pos[i][j] = (int **) malloc(size * sizeof(int *));
            for (int k = 0; k < size; k++) {
                pos[i][j][k] = (int *) malloc(3 * sizeof(int));
            }
        }
    }
    parts = (int **) malloc(npart * sizeof(int *));
    for (int i = 0; i < npart; i++) {
        parts[i] = (int *) malloc(6 * sizeof(int));
    }

    // Memory allocation for cluster analysis
    tempCluster = (int **) malloc(npart * sizeof(int *));
    for (int i = 0; i < npart; i++) {
        tempCluster[i] = (int *) malloc(6 * sizeof(int));
    }

    clusterList = (int ***) malloc(npart * sizeof(int **));
    for (int i = 0; i < npart; i++) {
        clusterList[i] = (int **) malloc(npart * sizeof(int *));
        for (int j = 0; j < npart; j++) {
            clusterList[i][j] = (int *) malloc(6 * sizeof(int));
        }
    }

    // Keeps a lookout for termination signals
    signal(SIGINT, intHandler);
    signal(SIGTERM, intHandler);
    signal(SIGSTOP, intHandler);


    //                              Simulation begins here...
    // -------------------------------------------------------------------------------------
    resetPos();
    molInit();
    initializeClusterList();
    resetTempCluster();

    // Internal energy check
    double energy = getTotalEnergy();
    double enCheck;

    printf("Size: %d\tTotal Particles: %d\n", size, npart);
    printf("Num Linkers: %d\n", npartLink);
    printf("Num Central: %d\n", npartMod);
    printf("CovE: %f\tVdwE: %f\n", CovE, VdwE);
    printf("\n");

    fprintf(outputFP, "Size: %d\tTotal Particles: %d\n", size, npart);
    fprintf(outputFP, "Num Linkers: %d\n", npartLink);
    fprintf(outputFP, "Num Central: %d\n", npartMod);
    fprintf(outputFP, "CovE: %f\tVdwE: %f\n\n", CovE, VdwE);

    printf("Step\tActualE  \tCalcE    \tDifference\n");
    fprintf(outputFP, "Step\tActualE  \tCalcE    \tDifference\n");
    while (keepRunning) {
        for (int i = 0; i < numSweeps; i++) {
            double dE = 0.0;

            // Check for consistency in energy calculations
            if (i % (numSweeps / 10) == 0) {
                enCheck = getTotalEnergy();
                printf("%d\t%f\t%f\t%f\n", i, energy, enCheck, (energy - enCheck));
                fprintf(outputFP, "%d\t%f\t%f\t%f\n", i, energy, enCheck, (energy - enCheck));
                fflush(stdout);
            }

            // Iterate through parts twice so each particle has equal probability of
            // performing a rotation and a translation
            for (int j = 0; j < npart * 2; j++) {
                // Randomly pick rotation or translation
                switch ((int) gsl_rng_uniform_int(rng, 2)) {
                    case 0:
                        // Randomly pick a particle to rotate
                        dE += proposeRotation((int) gsl_rng_uniform_int(rng, (unsigned long)npart));
                        break;
                    case 1:
                        // Randomly pick a particle to move
                        dE += proposeMove((int) gsl_rng_uniform_int(rng, (unsigned long)npart));
                        break;
                    default:
                        printf("ERROR: You've managed to have the number generator go beyond it's range\n");
                        printf("Save seed: %d and report to GSL.\n", seed);
                        break;
                }
            }
            energy += dE;

            // Print to xyz file every specified number of frames
            if ((i % (numSweeps / NUM_FRAMES)) == 0) {
                printXYZ();
            }
        } // End MC sweeps

        keepRunning = 0;
    } // End while (keepRunning)
    clusterFP = fopen("clusters.txt", "w");
    printf("\nCluster Search\n---------------------\n");
    findClusters();
    analyzeClusters();
    fclose(outputFP);
    fclose(clusterFP);
    fclose(xyzFP);
    return 0;
}

/*
 * -------------------------------------------------------------------------------------
 *                                    OUTPUT UTILITIES
 * -------------------------------------------------------------------------------------
 */
int timeStep = 0;

void printXYZ() {
    // Build header for xyz file
    fprintf(xyzFP, "ITEM: TIMESTEP\n%d\n", timeStep++);
    fprintf(xyzFP, "ITEM: NUMBER OF ATOMS\n%d\n", (4*npartMod + 3*npartLink));
    fprintf(xyzFP, "ITEM: BOX BOUNDS pp pp pp\n");
    for (int i = 0; i < 3; i++){
        fprintf(xyzFP,"-1 %d\n", size);
        fprintf(xyzFP,"0 %d\n", size-1);
    }
    fprintf(xyzFP, "ITEM: ATOMS type x y z\n");

    // Print each particle
    for (int i = 0; i < npart; i++){
        double x = (double) parts[i][XCOORD];
        double y = (double) parts[i][YCOORD];
        double z = (double) parts[i][ZCOORD];
        switch (parts[i][TYPE]){
            case MODULATOR:
                switch (parts[i][ROTATION]){
                    case 0:
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x-0.49, y, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x+0.49, y, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y-0.49, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y+0.49, z);
                        break;
                    case 1:
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x-0.49, y, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x+0.49, y, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y, z-0.49);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y, z+0.49);
                        break;
                    case 2:
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y, z-0.49);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y, z+0.49);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y-0.49, z);
                        fprintf(xyzFP, "Al\t%f\t%f\t%f\n", x, y+0.49, z);
                        break;
                    default:
                        break;
                }
                break;
            case LINKER:
                switch (parts[i][ROTATION] / 2){ // Rotations 0 and 1 will be aligned with x-axis and so on...
                    case 0: // Aligned with x-axis
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x-0.49, y, z);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y, z);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x+0.49, y, z);
                        break;
                    case 1: // y-axis
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y-0.49, z);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y, z);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y+0.49, z);
                        break;
                    case 2: // z-axis
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y, z-0.49);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y, z);
                        fprintf(xyzFP, "P\t%f\t%f\t%f\n", x, y, z+0.49);
                        break;
                    default:
                        break;
                }
                break;
            default:
                printf("ERROR: All particles must be 1 of 2 types.\n");
                break;
        }
    }
}

void printArgumentError(int numArgs, int errorType) {
    if (errorType == 0) {
        printf("\nIncorrect arguments. Expected 6, 4, or 3, but parsed %d. Please try again.\n\n", numArgs - 1);
    }
    if (errorType == 1) {
        printf("\nNumber of Monte Carlo sweeps must be over 1000. Input number %d.\n\n", numArgs);
    }
    printf("Can take up to 6 arguments:\n");
    printf("'*' indicate fields required for passing 4 arguments.\'*-' for 3 arguments.\nAll optional fields have defaults if not provided.\n\n");
    printf("(1) *- Size\n");
    printf("(2) *- Number of linker molecules\n");
    printf("(3) *- Number of central molecules\n");
    printf("(4)    Covalent bond energy (default: %.2f)\n", defaultCovE);
    printf("(5)    Nonspecific interaction energy (default: %.2f)\n", defaultVdwE);
    printf("(6) *  Number of Monte Carlo sweeps (must be >1000) (default: %d)\n\n\n", defaultNumSweeps);
}
