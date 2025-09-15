#ifndef STRUCT_H_
#define STRUCT_H_ 

typedef struct 
{
    double *coordinates;
} Point; 

typedef struct
{
    Point *Points;
    double *coordinates;
    int numPoints;
} Centroid; 
 
 /*int clustring_to_K(Point* all_points, Centroid* centroids_lst, int K, int d, int iter, int num_points, double eps); */
 #endif
