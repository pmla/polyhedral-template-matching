#ifndef NORMALIZE_VERTICES_H
#define NORMALIZE_VERTICES_H

void subtract_barycentre(int num, double *points, double (*normalized)[3]);
double normalize_vertices(int num, double *points, double (*normalized)[3]);

#endif

