/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011-12  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * posest demo. The program accepts a text file with the camera intrinsics K and
 * another containing quintuples of matching point coordinates (i.e., X Y Z  x y
 * where x y is a point in the image and X Y Z its corresponding 3D preimage)
 * and estimates the pose R, t mapping 3D points to 2D ones.
 * The corresponding camera matrix maps 3D points M to image projections m as
 * m ~ K*R*M + t
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "posest.h"

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

/* reads the 3x3 intrinsic calibration matrix contained in a file */
static void readCalibParams(char *fname, double K[9])
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
    K[0]=K[1]=K[2]=K[3]=K[4]=K[5]=K[6]=K[7]=K[8]=0.0;
    return;
  }

  ungetc(ch, fp);

  for(i=0; i<3; i++){
    if(fscanf(fp, "%lf%lf%lf\n", K, K+1, K+2)!=3){
      fprintf(stderr, "cannot read three numbers from row %d in file %s, exiting\n", i+1, fname);
      exit(1);
    }
    K+=3;
  }

  fclose(fp);
}

#define INL_PCENT 0.75

int main(int argc, char *argv[])
{
double (*pts2D)[2], (*pts3D)[3];
register int i;
int npts, noutl, *outidx=NULL;
int cstfunc;
char *icalfile, *matchesfile;
double K[9], P[NUM_PPARAMS], rms, rmeds;

clock_t start_time, end_time;

  /* arguments parsing */
  if(argc!=3){
    fprintf(stderr, "Usage: %s <K> <matched points>\n", argv[0]);
    exit(1);
  }

  icalfile=argv[1];
  matchesfile=argv[2];

  readCalibParams(icalfile, K);
  npts=readMatchingPoints(matchesfile, &pts2D, &pts3D);

#if 0
  for(i=0, rms=0.0; i<npts; ++i)
    rms+=fabs(pts3D[i][2]);
  if(fabs(rms)<1e-12){
    fprintf(stderr, "%s: third coordinates of all 3D points are zero, presumably from a calibration grid\n", argv[0]);
  }
#endif

#if 0
  for(i=0; i<npts; ++i){
    printf("%g %g %g  %g %g\n", pts3D[i][0], pts3D[i][1], pts3D[i][2], pts2D[i][0], pts2D[i][1]);
  }
#endif

#ifdef NEED_OUTLIERS
  if((outidx=(int *)malloc(npts*sizeof(int)))==NULL){
    fprintf(stderr, "Memory allocation request failed in main()\n");
    exit(1);
  }
#endif /* NEED_OUTLIERS */

  fprintf(stdout, "Camera pose estimation using %d image matches\n", npts);

  //cstfunc=POSEST_REPR_ERR_NO_NLN_REFINE; // no NL refinement
  //cstfunc=POSEST_OBJSPC_ERR_LHM; // use LHM
  cstfunc=POSEST_REPR_ERR_NLN_REFINE; // use the reprojection error
  //cstfunc=POSEST_REPR_ERR_NLN_MLSL_REFINE; // use the reprojection error & MLSL

  start_time=clock();
  if(fabs(K[0])>1e-10){ // known focal
    double rt[NUM_RTPARAMS];

    posest(pts2D, pts3D, npts, INL_PCENT, K, rt, NUM_RTPARAMS, cstfunc, outidx, &noutl, 1);
    end_time=clock();
    fprintf(stdout, "\nEstimated motion ([rv t]) [%d outliers, %.2lf%%]\n", noutl, (double)(100.0*noutl)/npts);
    for(i=0; i<NUM_RTPARAMS; ++i){
      fprintf(stdout, "%.7g ", rt[i]);
    }
    fprintf(stdout, "\n");
    posest_PfromKRt(P, K, rt);
  }
  else{ // if supplied K has zero focal, estimate it along with r & t
    double rtf[NUM_RTFPARAMS];

    posest(pts2D, pts3D, npts, INL_PCENT, K, rtf, NUM_RTFPARAMS, cstfunc, outidx, &noutl, 1);
    end_time=clock();
    fprintf(stdout, "\nEstimated motion & focal length [%d outliers, %.2lf%%]\n", noutl, (double)(100.0*noutl)/npts);
    for(i=0; i<NUM_RTFPARAMS; ++i){
      fprintf(stdout, "%.7g ", rtf[i]);
    }
    fprintf(stdout, "\n");

    K[0]=K[4]=rtf[NUM_RTFPARAMS-1];
    posest_PfromKRt(P, K, rtf);
  }
  fprintf(stdout, "\nElapsed time: %.2lf seconds, %.2lf msecs\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC,
                        ((double) (end_time - start_time)) / (CLOCKS_PER_SEC/1000.0));


  /* print diagnostics regarding reprojection errors */
  posest_RMS_RMedS(pts2D, pts3D, npts, P, &rms, &rmeds);
  fprintf(stdout, "\nReprojection RMS and RMedS errors for input points: %g %g\n", rms, rmeds);
  fflush(stdout);

  { double q1, q2, q3;

    posest_quartiles(pts2D, pts3D, npts, P, &q1, &q2, &q3);
    fprintf(stdout, "\t25%%, 50%% and 75%% quartiles: %g %g %g\n", q1, q2, q3);
    fflush(stdout);
  }

#ifdef NEED_OUTLIERS
  fprintf(stdout, "Indices of the %d outlying pairs:\n", noutl);
  for(i=0; i<noutl; ++i){
    fprintf(stdout, "%d ", outidx[i]);
    if(i && !(i%30)) fputc('\n', stdout);
  }
  fputc('\n', stdout);

  if(outidx) free(outidx);
#endif /* NEED_OUTLIERS */

  free(pts2D);
  free(pts3D);

  return 0;
}
