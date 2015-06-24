/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Restart PLUTO from binary or HDF5 double precision data files.
 
  This file collects the necessary functions for restarting PLUTO 
  from a double precision binary or HDF5 file in the static grid
  version of the code.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* *********************************************************************  */
void Restart (Input *ini, int nrestart, int type, Grid *grid)
/*!
 * Read input binary / hdf5 data.
 *
 * \param [in] ini
 * \param [in] nrestart
 * \param [in] type specifies the output data type (type should be 
 *             either DBL_OUTPUT or DBL_H5_OUTPUT).
 * \param [in]      grid  pointer to an array of Grid structures
 *
 ***********************************************************************  */
{
  int     i,j,k;
  int     nv, single_file, origin, nlines=0;
  int     swap_endian=0;
  char    fname[512], fout[512], str[512];
  double  dbl;
  void   *Vpt;
  Output *output;
  FILE   *fbin;

/* ----------------------------------------------------------
    get the pointer to the output format specified by "type"
   ---------------------------------------------------------- */

  for (nv = 0; nv < MAX_OUTPUT_TYPES; nv++){
    output = ini->output + nv;
    if (output->type == type) break;
  }

/* -------------------------------------------------------
    compare the endianity of the restart file (by reading
    the corresponding entry in dbl.out or dbl.h5.out) 
    with that of the current architecture.
    Turn swap_endian to 1 if they're different.
   ------------------------------------------------------- */

  if (prank == 0){
    if (type == DBL_OUTPUT) {
      sprintf (fout,"%s/dbl.out",ini->output_dir);
      fbin = fopen (fout, "r");
    } else if (type == DBL_H5_OUTPUT) {
      sprintf (fout,"%s/dbl.h5.out",ini->output_dir);
      fbin = fopen (fout, "r");
    }
    if (fbin == NULL){
      print1 ("! Restart: cannot find dbl.out or dbl.h5.out\n");
      QUIT_PLUTO(1);
    }

    while (fgets(str, 512, fbin) != 0) nlines++;  /* -- count lines in dbl.out -- */
    rewind(fbin);
    if (nrestart > nlines-1){
      printf ("! Restart: position too large\n");
      QUIT_PLUTO(1);
    }
    origin = (nrestart >= 0 ? nrestart:(nlines+nrestart));
    for (nv = origin; nv--;   ) while ( fgetc(fbin) != '\n'){}
    fscanf(fbin, "%d  %lf  %lf  %d  %s  %s\n",&nv, &dbl, &dbl, &nv, str, str);
    if ( (!strcmp(str,"big")    &&  IsLittleEndian()) ||
         (!strcmp(str,"little") && !IsLittleEndian())) {
      swap_endian = 1;
      print1 ("> Restart: endianity is reversed\n");
    }
    fclose(fbin);
  }
  #ifdef PARALLEL
   MPI_Bcast (&swap_endian, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

/* ---------------------------------------------
    read restart.out and get RunTime structure
   --------------------------------------------- */

  RestartGet (ini, nrestart, type, swap_endian);
  if (type == DBL_H5_OUTPUT){
    #ifdef USE_HDF5
     ReadHDF5 (output, grid);
    #endif
    return;
  }

  print1 ("> restarting from file #%d (dbl)\n",output->nfile);
  single_file = strcmp(output->mode,"single_file") == 0;
  
/* -----------------------------------------------------------------
           For .dbl output, read data from disk
   ----------------------------------------------------------------- */

  if (single_file){ 
    int  sz;
    long long offset;

    sprintf (fname, "%s/data.%04d.dbl", output->dir, output->nfile);
    offset = 0;
    #ifndef PARALLEL
     fbin = OpenBinaryFile (fname, 0, "r");
    #endif
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;

      if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
        sz = SZ;
        Vpt = (void *)output->V[nv][0][0];
      } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
        sz  = SZ_stagx;
        Vpt = (void *)(output->V[nv][0][0]-1);
      } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
        sz = SZ_stagy;
        Vpt = (void *)output->V[nv][0][-1];
      } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
         sz = SZ_stagz;
         Vpt = (void *)output->V[nv][-1][0];
      }
      #ifdef PARALLEL
       fbin = OpenBinaryFile (fname, sz, "r");
       AL_Set_offset(sz, offset);
      #endif
      ReadBinaryArray (Vpt, sizeof(double), sz, fbin,
                       output->stag_var[nv], swap_endian);
      #ifdef PARALLEL
       offset = AL_Get_offset(sz);
       CloseBinaryFile(fbin, sz);
      #endif
    }
    #ifndef PARALLEL
     CloseBinaryFile(fbin, sz);
    #endif

  }else{

    int  sz;
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;
      sprintf (fname, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                       output->nfile, output->ext);

      if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
        sz = SZ;
        Vpt = (void *)output->V[nv][0][0];
      } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
        sz  = SZ_stagx;
        Vpt = (void *)(output->V[nv][0][0]-1);
      } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
        sz = SZ_stagy;
        Vpt = (void *)output->V[nv][0][-1];
      } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
         sz = SZ_stagz;
         Vpt = (void *)output->V[nv][-1][0];
      }
      fbin = OpenBinaryFile (fname, sz, "r");
      ReadBinaryArray (Vpt, sizeof(double), sz, fbin,
                       output->stag_var[nv], swap_endian);
      CloseBinaryFile (fbin, sz);
    }
  }
}

static int counter = -1;

/* ********************************************************************* */
void RestartGet (Input *ini, int nrestart, int out_type,
                  int swap_endian)
/*!
 * Collect runtime information needed for (potential)
 * later restarts.
 *
 *********************************************************************** */
{
  int  origin, n, k;
  char fout[512];
  Runtime runtime;
  FILE *fr;

  if (nrestart < 0){
    printf ("! negative restart file temporarily disabled\n");
    QUIT_PLUTO(1);
  }

/* -------------------------------------------------
    Open "restart.out" and scan line by line until
    the output type specified by "out_type" has
    nfile = nrestart.
    counter will contain the line number where this
    occurs. Processor 0 does the actual reading.
   ------------------------------------------------- */

  if (prank == 0) {
    sprintf (fout,"%s/restart.out",ini->output_dir);
    fr = fopen (fout, "rb");
    if (fr == NULL){
      print1 ("! RestartGet: cannot find restart.out\n");
      QUIT_PLUTO(1);
    }

    origin = (nrestart < 0 ? SEEK_END:SEEK_SET);
    k = 0;
    while (counter == -1){
      if (feof(fr)){
        print("! RestartGet: end of file encountered.\n");
        QUIT_PLUTO(1);
      }
      fseek (fr, k*sizeof(Runtime), origin);
      fread (&runtime, sizeof (Runtime), 1, fr);
      for (n = 0; n < MAX_OUTPUT_TYPES; n++){
        if (swap_endian) SWAP_VAR(runtime.nfile[n]);
        if (ini->output[n].type == out_type && runtime.nfile[n] == nrestart){
          counter = k;
        }
      }
      k++;
    }
    fclose(fr);
    if (swap_endian){
      SWAP_VAR(runtime.t);
      SWAP_VAR(runtime.dt);
      SWAP_VAR(runtime.nstep);
    }
  }

/* printf ("counter = %d\n",counter); */

  #ifdef PARALLEL
   MPI_Bcast (&runtime, sizeof (Runtime), MPI_BYTE, 0, MPI_COMM_WORLD);
  #endif

  g_time = runtime.t;
  g_dt   = runtime.dt;
  g_stepNumber     = runtime.nstep;

/*printf ("Getting Runtime Structure\n"); */
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    ini->output[n].nfile = runtime.nfile[n];
/* if (runtime.nfile[n] >= 0)  printf ("output = %d, nfile = %d\n",n,runtime.nfile[n]);  */
  }
}
/* ********************************************************************* */
void RestartDump (Input *ini)
/*!
 * Write runtime information needed for (potential) later restarts.
 *
 *********************************************************************** */
{
  int n;
  char fout[512];
  Runtime runtime;
  FILE *fr;

/* --------------------------------------------------
    define runtime structure elements here
   -------------------------------------------------- */

  runtime.t  = g_time;
  runtime.dt = g_dt;
  runtime.nstep = g_stepNumber;
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    runtime.nfile[n] = ini->output[n].nfile;
/* if (runtime.nfile[n] >= 0) printf ("output = %d, nfile = %d\n",n,runtime.nfile[n]); */
  }  

/* --------------------------------------------------
         dump structure to disk
   -------------------------------------------------- */

  counter++;
/* printf ("Dumping Runtime Structure; counter = %d\n", counter);  */
  if (prank == 0) {
    sprintf (fout,"%s/restart.out",ini->output_dir);
    if (counter == 0) {
      fr = fopen (fout, "wb");
    }else {
      fr = fopen (fout, "r+b");
      fseek (fr, counter*sizeof(Runtime), SEEK_SET); 
    }

    fwrite (&runtime, sizeof(Runtime), 1, fr);
    fclose(fr);
  }
}
