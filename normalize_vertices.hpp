#ifndef NORMALIZE_VERTICES_HPP
#define NORMALIZE_VERTICES_HPP

void subtract_barycentre(int num, double (*points)[3], double (*normalized)[3]);
double normalize_vertices(int num, double (*points)[3], double (*normalized)[3]);

#endif

