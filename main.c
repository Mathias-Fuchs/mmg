#include <stdio.h>
#include <stdlib.h>
#include <mmg/mmg2d/libmmg2d.h>
#include <math.h>

int main() {
  MMG5_pMesh mesh2 = NULL;
  MMG5_pSol met2 = NULL;
  int ier = 0;
  MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh2, MMG5_ARG_ppMet, &met2, MMG5_ARG_end);
  MMG2D_Set_iparameter(mesh2, met2, MMG2D_IPARAM_verbose, 10);

  // 4 vertices, no tris, no quads, no edges so far
  MMG2D_Set_meshSize(mesh2, 4, 0, 0, 0);

  // a square [0, 2 * M_PI] x [-M_PI / 2, M_PI / 2] for spherical coordinates in usual lon/lat way
  double v[8] = {0, -M_PI/2, 2 * M_PI, -M_PI/2, 2 * M_PI, M_PI/2, 0, M_PI/2};
  MMG2D_Set_vertices(mesh2, v, NULL);
  
  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hsiz, 0.01);

  // generate a regular fine mesh of the unit square in meshing mode
  ier = MMG2D_mmg2dmesh(mesh2, met2);
  if (ier) {
    fprintf(stderr, "error %i during meshing.\n", ier);
    exit(1);
  }

  // save the "computational geometry" mesh
  MMG2D_saveMshMesh(mesh2, NULL, "cg.msh");

  
  // remesh with anisotropic metric
  int np, nt, nquad, na;
  MMG2D_Get_meshSize(mesh2, &np, &nt, &nquad, &na);
  double* verts = (double*) malloc(2 * np * sizeof(double));
  MMG2D_Get_vertices(mesh2, verts, NULL, NULL, NULL);
  MMG2D_Set_solSize(mesh2, met2, MMG5_Vertex, np, MMG5_Tensor);
  for (int i = 0; i < np; i++) {

    // latitude
    double y = verts[2 * i + 1];
    // metric on the sphere, see standard textbooks on elementary differential geometry
    // Gaussian fundamental quantities
    double E = cos(y) * cos(y);
    double F = 0.0;
    double G = 1.0;
    
    // the following factor was found by trial and error
    // why does it not work with a factor of one?
    // how to compute the factor in real life?
    double factor = 1000.0;

    // we add a small amount to the 1-1-entry to make the metric non-singular everywhere
    MMG2D_Set_tensorSol(met2, factor * E + 0.1, factor * F, factor * G, i+1);
  }
  
  // disable hsiz because it is incompatible with an anisotropic metric
  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hsiz, -1);

  // set gradation to a not too restrictive value
  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hgrad, 1.5);

  // do it, remesh!
  ier = MMG2D_mmg2dlib(mesh2, met2);

  if (ier != 0) {
    fprintf(stdout, "error %i \n", ier);
    exit(1);
  }
  
  MMG2D_saveMshMesh(mesh2, NULL, "out2.msh");

  // map to 3d and save as obj file
  MMG2D_Get_meshSize(mesh2, &np, &nt, &nquad, &na);
  verts = realloc(verts, 2 * np * sizeof(double));
  int* tris = (int*) malloc(3 * nt * sizeof(int));
  MMG2D_Get_vertices(mesh2, verts, NULL, NULL, NULL);
  MMG2D_Get_triangles(mesh2, tris, NULL, NULL);
  
  FILE *fp = fopen("sphere.obj", "w+");
  for (int i = 0; i < np; i++) {
    double x = verts[2 * i];
    double y = verts[2 * i + 1];
    double newx = cos(x) * cos(y);
    double newy = sin(x) * cos(y);
    double newz = sin(y);
    fprintf(fp, "v %f %f %f\n", newx, newy, newz);
  }
  for (int i = 0; i < nt; i++)
    fprintf(fp, "f %i %i %i\n", tris[3 * i], tris[3 * i +1], tris[3 * i + 2]);

  fclose(fp);
  free(verts);
  free(tris);
  return 0;
}
