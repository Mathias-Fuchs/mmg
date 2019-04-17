/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file common/mmg.c
 * \brief Common part for functions used in mmgs.c and mmg3d.c files.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 04 2015
 * \copyright GNU Lesser General Public License.
 **/

#include "mmgcommon.h"

/**
 * \param *prog pointer toward the program name.
 *
 * Print help for common options of mmg3d and mmgs.
 *
 */
void MMG5_mmgUsage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h        Print this message\n");
  fprintf(stdout,"-v [n]    Tune level of verbosity, [-1..10]\n");
  fprintf(stdout,"-m [n]    Set maximal memory size to n Mbytes\n");
  fprintf(stdout,"-d        Turn on debug mode\n");
  fprintf(stdout,"-val      Print the default parameters values\n");
  fprintf(stdout,"-default  Save a local parameters file for default parameters"
          " values\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution or metric file\n");

  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-ar     val  angle detection\n");
  fprintf(stdout,"-nr          no angle detection\n");
  fprintf(stdout,"-hmin   val  minimal mesh size\n");
  fprintf(stdout,"-hmax   val  maximal mesh size\n");
  fprintf(stdout,"-hsiz   val  constant mesh size\n");
  fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
  fprintf(stdout,"-hgrad  val  control gradation\n");
  fprintf(stdout,"-ls     val  create mesh of isovalue val (0 if no argument provided)\n");

}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG5_mmgDefaultValues(MMG5_pMesh mesh) {

  fprintf(stdout,"\nDefault parameters values:\n");

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"verbosity                 (-v)      : %d\n",
          mesh->info.imprim);
  fprintf(stdout,"maximal memory size       (-m)      : %zu MB\n",
          mesh->memMax/MMG5_MILLION);


  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"angle detection           (-ar)     : %lf\n",
          180/M_PI*acos(mesh->info.dhd) );
  fprintf(stdout,"minimal mesh size         (-hmin)   : 0.001 of "
          "the mesh bounding box if no metric is provided, 0.1 times the "
          "minimum of the metric sizes otherwise.\n");
  fprintf(stdout,"maximal mesh size         (-hmax)   : size of "
          "the mesh bounding box without metric, 10 times the maximum of the "
          "metric sizes otherwise.\n");
  fprintf(stdout,"Hausdorff distance        (-hausd)  : %lf\n",
          mesh->info.hausd);
  fprintf(stdout,"gradation control         (-hgrad)  : %lf\n",
          exp(mesh->info.hgrad));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \return npar, the number of local parameters at triangles if success,
 * 0 otherwise.
 *
 * Count the local default values at triangles and fill the list of the boundary
 * references.
 *
 */
inline
int MMG5_countLocalParamAtTri( MMG5_pMesh mesh,MMG5_iNode **bdryRefs) {
  int         npar,k,ier;

  /** Count the number of different boundary references and list it */
  (*bdryRefs) = NULL;

  k = mesh->nt? mesh->tria[1].ref : 0;

  /* Try to alloc the first node */
  ier = MMG5_Add_inode( mesh, bdryRefs, k );
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate the first boundary"
           " reference node.\n",__func__);
    return 0;
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->nt; ++k ) {
    ier = MMG5_Add_inode( mesh, bdryRefs, mesh->tria[k].ref );

    if ( ier < 0 ) {
      printf("  ## Warning: %s: unable to list the tria references."
             " Uncomplete parameters file.\n",__func__ );
      break;
    }
    else if ( ier ) ++npar;
  }

  return npar;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \param npar number of local param at triangles.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write the local default values at triangles in the parameter file.
 *
 */
inline
int MMG5_writeLocalParamAtTri( MMG5_pMesh mesh, MMG5_iNode *bdryRefs,
                                FILE *out ) {
  MMG5_iNode *cur;

  cur = bdryRefs;
  while( cur ) {
    fprintf(out,"%d Triangle %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  MMG5_Free_ilinkedList(mesh,bdryRefs);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param mesh pointer toward the msh value.
 *
 * Update the msh value if we detect that the user want to force output at Gmsh
 * or Medit format.
 *
 */
void MMG5_chooseOutputFormat(MMG5_pMesh mesh, int *msh) {
  int len;

  len = strlen(mesh->nameout);

  if ( ( len>4 && !strcmp(&mesh->nameout[len-5],".mesh") ) ||
       ( len>5 && !strcmp(&mesh->nameout[len-6],".meshb") ) )
    *msh = 0;
  else if ( ( len>3 && !strcmp(&mesh->nameout[len-4],".msh") ) ||
            ( len>4 && !strcmp(&mesh->nameout[len-5],".mshb") ))
    *msh = 1;
  else
    *msh = 0;

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 */
void MMG5_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      isqhmin, isqhmax;
  int         i,k,iadr,sethmin,sethmax;

  assert ( mesh->info.optim || mesh->info.hsiz > 0. );

  /* Detect the point used only by prisms */
  if ( mesh->nprism ) {
    for (k=1; k<=mesh->np; k++) {
      mesh->point[k].flag = 1;
    }
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      for (i=0; i<4; i++) {
        mesh->point[pt->v[i]].flag = 0;
      }
    }
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  sethmin = sethmax = 1;
  if ( mesh->info.hmin < 0 ) {
    sethmin = 0;
    if ( met->size == 1 ) {
      mesh->info.hmin = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmin = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+3]);
        mesh->info.hmin = MG_MAX(mesh->info.hmin,met->m[iadr+5]);
      }
      mesh->info.hmin = 1./sqrt(mesh->info.hmin);
    }
  }

  if ( mesh->info.hmax < 0 ) {
    sethmax = 0;
    if ( met->size == 1 ) {
      mesh->info.hmax = 0.;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
      }
    }
    else if ( met->size == 6 ){
      mesh->info.hmax = FLT_MAX;
      for (k=1; k<=mesh->np; k++)  {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || ppt->flag ) continue;
        iadr = met->size*k;
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+3]);
        mesh->info.hmax = MG_MIN(mesh->info.hmax,met->m[iadr+5]);
      }
      mesh->info.hmax = 1./sqrt(mesh->info.hmax);
    }
  }


  if ( !sethmin ) {
    mesh->info.hmin *=.1;
    /* Check that user has not given a hmax value lower that the founded
     * hmin. */
    if ( mesh->info.hmin > mesh->info.hmax ) {
      mesh->info.hmin = 0.1*mesh->info.hmax;
    }
  }
  if ( !sethmax ) {
    mesh->info.hmax *=10.;
    /* Check that user has not given a hmin value bigger that the founded
     * hmax. */
    if ( mesh->info.hmax < mesh->info.hmin ) {
      mesh->info.hmax = 10.*mesh->info.hmin;
    }
  }

  /* vertex size */
  if ( met->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
    }
  }
  else if ( met->size == 6 ) {
    isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
    isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      iadr = 6*k;
      met->m[iadr]   = MG_MAX(isqhmax,MG_MIN(isqhmin,met->m[iadr]));
      met->m[iadr+3] = met->m[iadr];
      met->m[iadr+5] = met->m[iadr];
    }
  }

  return;
}
