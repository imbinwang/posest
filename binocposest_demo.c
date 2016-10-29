/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011-13  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * Binocular posest demo. The program accepts text files with the two camera
 * matrices and matching point coordinates in each image, in the posest_demo
 * format (i.e., X Y Z  x y  for 3D points and their 2D projections)
 * It estimates the pose R, t mapping 3D points to 2D ones. The pose is computed
 * having the left (first) camera as origin.
 * The corresponding camera matrix maps 3D points M to image projections m as
 * mL ~ KL*R*M + t
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "sam.h"
#include "compiler.h"
#include "posest.h"
#include "poseproj.h"


#undef NEED_OUTLIERS // define this to print information about the detected outliers

#define MAXSTRLEN 1024

/* read matching points from a file */
static int readMatchingPoints(char *fname, double (**pts2D)[2], double (**pts3D)[3])
{
register int i;
int ncoords, nmatches;
double coords[5];
FILE *fp;
char buf[MAXSTRLEN], *ptr;
long int fpos;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s\n", fname);
    exit(1);
  }

  do{
    fpos=ftell(fp);
    ptr=fgets(buf, MAXSTRLEN-1, fp);
    if(!ptr || ferror(fp)){
      fprintf(stderr, "File %s: error reading line \"%s\"\n", fname, ptr);
      exit(1);
    }
  } while(!feof(fp) && buf[0]=='#'); /* skip comments */

  if(feof(fp)){
    fclose(fp);
    return 0;
  }

  ncoords=sscanf(buf, "%lf%lf%lf%lf%lf", coords, coords+1, coords+2, coords+3, coords+4);
  if(ncoords==5){ /* no lines number */
    for(nmatches=1; !feof(fp); nmatches++){
      i=fscanf(fp, "%*g%*g%*g%*g%*g\n");
      if(ferror(fp)){
        fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, nmatches + 1);
        exit(1);
      }
    }

    if(fseek(fp, fpos, SEEK_SET)){ /* rewind right after any comment lines */
      fprintf(stderr, "fseek failed in readMatchingPoints()\n");
      exit(1);
    }
  }
  else{
    sscanf(buf, "%d", &nmatches);
  }

  *pts2D=(double (*)[2])malloc(nmatches*sizeof(double[2]));
  *pts3D=(double (*)[3])malloc(nmatches*sizeof(double[3]));
  if(!pts2D || !pts3D){
    fprintf(stderr, "Memory allocation request failed in readMatchingPoints()\n");
    exit(1);
  }

  /* read in points and store them */
  for(i=0; !feof(fp); i++){
    ncoords=fscanf(fp, "%lf%lf%lf%lf%lf\n", (*pts3D)[i], (*pts3D)[i]+1, (*pts3D)[i]+2, (*pts2D)[i], (*pts2D)[i]+1);
    if(ncoords==EOF) break;

    if(ncoords!=5){
      fprintf(stderr, "File %s: line %d contains only %d coordinates\n", fname, i + 1, ncoords);
      exit(1);
    }

    if(ferror(fp)){
      fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, i + 1);
      exit(1);
    }

  }
  fclose(fp);

  if(i!=nmatches){
    fprintf(stderr, "number of actuall points in file %s does not agree with that in first line (%d != %d)!\n",
                     fname, i, nmatches);
    exit(1);
  }

  return nmatches;
}

/* reads the 3x4 camera matrix contained in a file */
static void readPMat(char *fname, double P[12])
{
  FILE *fp;
  int i, ch=EOF;
  char buf[MAXSTRLEN], *ptr;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", fname);
    exit(1);
  }

  while(!feof(fp) && (ch=fgetc(fp))=='#') /* skip comments */
    ptr=fgets(buf, MAXSTRLEN-1, fp);

  if(feof(fp)){
    fclose(fp);
    P[0]=P[1]=P[2]=P[3]=P[4]=P[5]=
    P[6]=P[7]=P[8]=P[9]=P[10]=P[11]=0.0;
    return;
  }

  ungetc(ch, fp);

  for(i=0; i<3; i++){
    if(fscanf(fp, "%lf%lf%lf%lf\n", P, P+1, P+2, P+3)!=4){
      fprintf(stderr, "cannot read four numbers from row %d in file %s, exiting\n", i+1, fname);
      exit(1);
    }
    P+=4;
  }

  fclose(fp);
}


#define INLPCENT  0.4

int main(int argc, char *argv[])
{
int nptsL;
double (*pts2DL)[2], (*pts3DL)[3];
double PextL[NUM_PPARAMS];

int nptsR;
double (*pts2DR)[2], (*pts3DR)[3];
double PextR[NUM_PPARAMS];

double rts[NUM_RTPARAMS+1]; // r, t, scale
double inscale=1.0;
int howto, noutl, *outidx=NULL;
int verbose=1, doscale=0;
//double KL[9]={3314.14, 3.63798e-012, 675.401, 0, 3328.59, 474.463, 0, 0, 1};
//double KR[9]={3320.97, 0, 698.217, 0, 3329.69, 458.488, 0, 0, 1};
double q1, q2, q3;
clock_t start_time, end_time;

  /* arguments parsing */
  if(argc!=5){
    fprintf(stderr, "Usage: %s <Pl> <left points> <Pr> <right points>\n", argv[0]);
    exit(1);
  }

  readPMat(argv[1], PextL);
  nptsL=readMatchingPoints(argv[2], &pts2DL, &pts3DL);

  readPMat(argv[3], PextR);
  nptsR=readMatchingPoints(argv[4], &pts2DR, &pts3DR);

#ifdef NEED_OUTLIERS
  if((outidx=(int *)malloc((nptsL+nptsR)*sizeof(int)))==NULL){
    fprintf(stderr, "Memory allocation request failed in main()\n");
    exit(1);
  }
#endif /* NEED_OUTLIERS */

//nmatches=nptsL; posest(pts2DL, pts3DL, nptsL, INLPCENT, KL, rts, NUM_RTPARAMS, POSEST_REPR_ERR_NLN_REFINE, NULL, &noutl, verbose); // CHECKME

  //howto=POSEST_REPR_ERR_NLN_REFINE;
  howto=POSEST_REPR_ERR_NLN_MLSL_REFINE;

  start_time=clock();
  posestBinoc(pts2DL, pts3DL, nptsL, PextL, pts2DR, pts3DR, nptsR, PextR,
              INLPCENT, inscale, rts, doscale, howto, outidx, &noutl, verbose);
  end_time=clock();
  printf("Refined motion ([rv t  scale]) [%d outliers, %.2lf%%]\n", noutl, (double)(100.0*noutl)/(nptsL+nptsR));
  printf("\t%g %g %g %g %g %g  %.6lf\n", rts[0], rts[1], rts[2], rts[3], rts[4], rts[5], rts[6]);
  printf("\n");
  fprintf(stdout, "\nElapsed time: %.2lf seconds, %.2lf msecs\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC,
                              ((double) (end_time - start_time)) / (CLOCKS_PER_SEC/1000.0));

  //posestBinoc_RMS_RMedS(pts2DL, pts3DL, nptsL, PextL, pts2DR, pts3DR, nptsR, PextR, rts, &rms, &rmeds);
  //fprintf(stdout, "\nReprojection RMS and RMedS errors for input points: %g %g\n", rms, rmeds);

  posestBinoc_quartiles(pts2DL, pts3DL, nptsL, PextL, pts2DR, pts3DR, nptsR ,PextR, rts, &q1, &q2, &q3);
  fprintf(stdout, "Reprojection errors 25%%, 50%% and 75%% quartiles: %g %g %g\n", q1, q2, q3);

#ifdef NEED_OUTLIERS
  {
  int i;

  fprintf(stdout, "Indices of the %d outlying pairs:\n", noutl);
  for(i=0; i<noutl; ++i){
    fprintf(stdout, "%d ", outidx[i]);
    if(i && !(i%30)) fputc('\n', stdout);

    //if(outidx[i]<nptsL) then outidx[i] is an outlier from the left points
    //else outidx[i]-nptsL is an outlier from the right points
  }
  fputc('\n', stdout);

  if(outidx) free(outidx);
  }
#endif /* NEED_OUTLIERS */

  free(pts2DL); free(pts3DL);
  free(pts2DR); free(pts3DR);

  return 0;
}
