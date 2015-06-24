/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Read runtime information from pluto.ini.

  Parse and read runtime information from the initialization file 
  pluto.ini (default) and sets value of the input structure.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    May 29, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void GetOutputFrequency(Output *, char *);

#define NOPT      32              /*  # of possible options in a menu */
#define NLEN      128              /*  # default string length         */

#define COMPARE(s1,s2,ii) \
        for ( (ii) = 1 ; (ii) < NOPT && !(strcmp ( (s1), (s2)) == 0); (ii)++);

/* ********************************************************************* */
int Setup (Input *input, Cmd_Line *cmd_line, char *ini_file)
/*!
 * Open and parse the initialization file. 
 * Assign values to the input structure.
 *
 * \param [out]  input     pointer to an Input structure
 * \param [in]   cmd_line  pointer to a Cmd_Line structure (useful, e.g.,
 *                         to resize the domain using the \c -xres option)
 * \param [in]   ini_file  the name of the initialization file (default
 *                         is "pluto.ini") specified with the \c -i option.
 *
 *********************************************************************** */
{
  int    idim, ip, ipos, itype, nlines;
  char   *bound_opt[NOPT], str_var[512], *str;
  char  *glabel[]     = {"X1-grid", "X2-grid","X3-grid"};
  char  *bbeg_label[] = {"X1-beg", "X2-beg","X3-beg"};
  char  *bend_label[] = {"X1-end", "X2-end","X3-end"};
  double dbl_var, rx;
  Output *output;
  FILE *fp;

  for (itype = 0; itype < NOPT; itype++) {
    bound_opt[itype] = "0000";
  }

/*  ---------------------------------------------------
      available options are given as two set of names;
      This facilitates when updating the code and
      people are too lazy to read the manual !
    --------------------------------------------------- */

  bound_opt[OUTFLOW]      = "outflow";
  bound_opt[REFLECTIVE]   = "reflective";
  bound_opt[AXISYMMETRIC] = "axisymmetric";
  bound_opt[EQTSYMMETRIC] = "eqtsymmetric";
  bound_opt[PERIODIC]     = "periodic";
  bound_opt[SHEARING]     = "shearingbox";
  bound_opt[USERDEF]      = "userdef";

  input->log_freq = 1; /* -- default -- */
 
  nlines = ParamFileRead(ini_file);

/* ------------------------------------------------------------
                        [Grid] Section 
   ------------------------------------------------------------ */

  for (idim = 0; idim < 3; idim++){
    input->npatch[idim] = atoi(ParamFileGet(glabel[idim], 1));
    input->npoint[idim] = 0;

    ipos = 1;
    for (ip = 1; ip <= input->npatch[idim]; ip++) {

      input->patch_left_node[idim][ip] = atof(ParamFileGet(glabel[idim], ++ipos));
      input->patch_npoint[idim][ip]    = atoi(ParamFileGet(glabel[idim], ++ipos));
      input->npoint[idim]             += input->patch_npoint[idim][ip];
      input->grid_is_uniform[idim]     = 0;

      strcpy (str_var, ParamFileGet(glabel[idim], ++ipos)); 
/*
printf ("%f  %d %s\n",input->patch_left_node[idim][ip],input->patch_npoint[idim][ip],str_var);
*/
      if (strcmp(str_var,"u") == 0 || strcmp(str_var,"uniform") == 0) {
        input->patch_type[idim][ip] = UNIFORM_GRID;
        if (input->npatch[idim] == 1) input->grid_is_uniform[idim] = 1;        
      }else if (strcmp(str_var,"s") == 0 || strcmp(str_var,"strecthed") == 0) { 
        input->patch_type[idim][ip] = STRETCHED_GRID;
      }else if (strcmp(str_var,"l+") == 0){
        input->patch_type[idim][ip] = LOGARITHMIC_INC_GRID;
      }else if (strcmp(str_var,"l-") == 0){
        input->patch_type[idim][ip] = LOGARITHMIC_DEC_GRID;
      }else{ 
        printf ("\nSetup: You must specify either 'u', 's', 'l+' or 'l-' as grid-type in %s\n",
                ini_file);
        QUIT_PLUTO(1);
      }
    }
    
    input->patch_left_node[idim][ip] = atof(ParamFileGet(glabel[idim], ++ipos));

    if ( (ipos+1) != (input->npatch[idim]*3 + 3)) {
      printf ("! Setup: domain #%d setup is not properly defined \n", idim);
      QUIT_PLUTO(1);
    }
    if (idim >= DIMENSIONS && input->npoint[idim] != 1) {
      printf ("! Setup: %d point(s) on dim. %d is NOT valid, resetting to 1\n",
              input->npoint[idim],idim+1);
      input->npoint[idim]          = 1;
      input->npatch[idim]          = 1;
      input->patch_npoint[idim][1] = 1;
    }
  }

/* ------------------------------------------------------------
      Change the resolution if cmd_line->xres has been given
   ------------------------------------------------------------ */

  if (cmd_line->xres > 1) {
    rx =  (double)cmd_line->xres/(double)input->patch_npoint[IDIR][1];
    for (idim = 0; idim < DIMENSIONS; idim++){
      if (input->npatch[idim] > 1){  
        printf ("! Setup: -xres option works on uniform, single patch grid\n");
        QUIT_PLUTO(1);
      }
      
      dbl_var = (double)input->patch_npoint[idim][1];
      input->patch_npoint[idim][1] = MAX( (int)(dbl_var*rx), 1);
      dbl_var = (double)input->npoint[idim];
      input->npoint[idim] = MAX( (int)(dbl_var*rx), 1); 
    }  
  }

/* ------------------------------------------------------------
                     [Time] Section 
   ------------------------------------------------------------ */

  input->cfl         = atof(ParamFileGet("CFL", 1));

  if (ParamExist ("CFL_par")) input->cfl_par = atof(ParamFileGet("CFL_par", 1));
  else                            input->cfl_par = 0.8/(double)DIMENSIONS;

  if (ParamExist ("rmax_par")) input->rmax_par = atof(ParamFileGet("rmax_par", 1));
  else                             input->rmax_par = 100.0;

  input->cfl_max_var = atof(ParamFileGet("CFL_max_var", 1));
  input->tstop       = atof(ParamFileGet("tstop", 1));
  input->first_dt    = atof(ParamFileGet("first_dt", 1));

/* ------------------------------------------------------------
                     [Solver] Section 
   ------------------------------------------------------------ */

  sprintf (input->solv_type,"%s",ParamFileGet("Solver",1));

/* ------------------------------------------------------------
                     [Boundary] Section 
   ------------------------------------------------------------ */

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bbeg_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! Setup: don't know how to put left boundary '%s'  \n", str);
      QUIT_PLUTO(1);
    }
    input->lft_bound_side[idim] = itype;
  }

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bend_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! Setup: don't know how to put left boundary '%s'  \n", str);
      QUIT_PLUTO(1);
    }
    input->rgt_bound_side[idim] = itype;
  }

/* ------------------------------------------------------------
                     [Output] Section 
   ------------------------------------------------------------ */

  input->user_var = atoi(ParamFileGet("uservar", 1));
  for (ip = 0; ip < input->user_var; ip++){

    if ( (str = ParamFileGet("uservar", 2 + ip)) != NULL){
      sprintf (input->user_var_name[ip], "%s", str);
    }else{
      printf ("! Setup: missing name after user var name '%s'\n", 
              input->user_var_name[ip-1]);
      QUIT_PLUTO(1);
    } 
  }

/* ---- set output directory ---- */

  sprintf (input->output_dir, "%s","./");  /* default value is current directory */
  if (ParamExist("output_dir")){
    str = ParamFileGet("output_dir",1);
    sprintf (input->output_dir, "%s",str);
  }
  SetOutputDir(input->output_dir);

/* -- check if we have write access and if the directory exists -- */

  sprintf (str_var,"%s/tmp09123.txt",input->output_dir); /* -- test file -- */
  fp = fopen(str_var,"w");  /* -- open test file for writing -- */
  if (fp == NULL){
    printf ("! Setup: cannot access directory '%s'.\n", input->output_dir);
    printf ("!        Please check that the directory exists\n");
    printf ("!        and you have write permission.\n");
    QUIT_PLUTO(1);
  }else{
    fclose(fp);
    remove(str_var);  /* -- remove file -- */
  }

/* ---- dbl output ---- */

  ipos = 0;
  output = input->output + (ipos++);
  output->type  = DBL_OUTPUT;
  output->cgs   = 0;  /* cannot write .dbl using cgs units */
  GetOutputFrequency(output, "dbl");

  sprintf (output->mode,"%s",ParamFileGet("dbl",3));
  #ifdef USE_ASYNC_IO
   if (    strcmp(output->mode,"single_file") 
        && strcmp(output->mode,"single_file_async")
        && strcmp(output->mode,"multiple_files")){
      printf ("! Setup: expecting 'single_file', 'single_file_async' ");
      printf ("or 'multiple_files' in dbl output\n");
      QUIT_PLUTO(1);
   }
  #else
   if (   strcmp(output->mode,"single_file")
       && strcmp(output->mode,"multiple_files")){
      printf (
      "! Setup: expecting 'single_file' or 'multiple_files' in dbl output\n");
      QUIT_PLUTO(1);
   }     
  #endif

 /* ---- flt output ---- */

  if (ParamExist("flt")){
    output = input->output + (ipos++);
    output->type  = FLT_OUTPUT;
    GetOutputFrequency(output, "flt");

    sprintf (output->mode,"%s",ParamFileGet("flt",3));  
    #ifdef USE_ASYNC_IO
     if (    strcmp(output->mode,"single_file") 
          && strcmp(output->mode,"single_file_async")
          && strcmp(output->mode,"multiple_files")){
        printf ("! Setup: expecting 'single_file', 'single_file_async' ");
        printf ("or 'multiple_files' in flt output\n");
        QUIT_PLUTO(1);
     }
    #else
     if (    strcmp(output->mode,"single_file") 
          && strcmp(output->mode,"multiple_files")){
        printf (
        "! Setup: expecting 'single_file' or 'multiple_files' in flt output\n");
        QUIT_PLUTO(1);
     }  
    #endif
    if (ParamFileHasBoth ("flt","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- hdf5 output -- */

  if (ParamExist("dbl.h5")){
    output = input->output + (ipos++);
    output->type  = DBL_H5_OUTPUT;
    output->cgs   = 0;  /* cannot write .h5 using cgs units */
    GetOutputFrequency(output, "dbl.h5");
  }
  if (ParamExist("flt.h5")){
    output = input->output + (ipos++);
    output->type  = FLT_H5_OUTPUT;
    output->cgs   = 0;  /* cannot write .h5 using cgs units */
    GetOutputFrequency(output, "flt.h5");
  }

 /* -- vtk output -- */

  if (ParamExist ("vtk")){
    output = input->output + (ipos++);
    output->type  = VTK_OUTPUT;
    GetOutputFrequency(output, "vtk");

    if (ParamFileGet("vtk",3) == NULL){
      printf ("! Setup: extra field missing in vtk output\n");
      QUIT_PLUTO(1);
    }
    sprintf (output->mode,"%s",ParamFileGet("vtk",3));
    if (   strcmp(output->mode,"single_file")
        && strcmp(output->mode,"multiple_files")){
       printf ("! Setup: expecting 'single_file' or 'multiple_files' in\n");
       printf ("         vtk output\n");
       QUIT_PLUTO(1);
    }
    if (ParamFileHasBoth ("vtk","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- tab output -- */

  if (ParamExist ("tab")){
    output = input->output + (ipos++);
    output->type  = TAB_OUTPUT;
    GetOutputFrequency(output, "tab");
    if (ParamFileHasBoth ("tab","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- ppm output -- */

  if (ParamExist ("ppm")){
    output = input->output + (ipos++);
    output->type  = PPM_OUTPUT;
    output->cgs   = 0;   /* Cannot write ppm in cgs units */
    GetOutputFrequency(output, "ppm");
  }

 /* -- png output -- */

  if (ParamExist ("png")){
    output = input->output + (ipos++);
    output->type  = PNG_OUTPUT;
    output->cgs   = 0;   /* Cannot write png in cgs units */
    GetOutputFrequency(output, "png");
  }


 /* -- log frequency -- */

  input->log_freq = atoi(ParamFileGet("log", 1));
  input->log_freq = MAX(input->log_freq, 1);
  
 /* -- set default for remaining output type -- */

  while (ipos < MAX_OUTPUT_TYPES){
    output = input->output + ipos;
    output->type   = -1;
    output->cgs    =  0;
    output->dt     = -1.0;
    output->dn     = -1;
    output->dclock = -1.0;
    ipos++;
  }

 /* -- analysis -- */

  if (ParamExist ("analysis")){
    input->anl_dt = atof(ParamFileGet("analysis", 1));
    input->anl_dn = atoi(ParamFileGet("analysis", 2));
  }else{
    input->anl_dt = -1.0;   /* -- defaults -- */
    input->anl_dn = -1;
  }

/* ------------------------------------------------------------
                   [Parameters] Section 
   ------------------------------------------------------------ */

  fp = fopen(ini_file,"r");
  
/* -- find position at "[Parameters" -- */

  for (ipos = 0; ipos <= nlines; ipos++){ 
    fgets(str_var, 512, fp);
    
    if (strlen(str_var) > 0) {
      str = strtok (str_var,"]");
      if (strcmp(str,"[Parameters") == 0) break;
    }
  }

  fgets(str_var, 512, fp); 
  
  for (ip = 0; ip < USER_DEF_PARAMETERS; ip++){
    fscanf (fp,"%s \n", str_var);
    dbl_var = atof(ParamFileGet(str_var,1));
    input->aux[ip] = dbl_var;
    fgets(str_var, sizeof(str_var), fp); /* use fgets to advance to next line */
  }
  fclose(fp);

  return(0);
}
#undef COMPARE
#undef NOPT
#undef NLEN

/* ********************************************************************* */
void GetOutputFrequency(Output *output, char *output_format)
/*!
 *  Set the intervals between output files. 
 *  This can be done in three different ways:
 *
 *  - dt: time interval in code units
 *  - dn: step interval
 *  - dclock: actual clock time (in hours)
 * 
 * However, dn and dclock are mutually exclusive.
 *
 *********************************************************************** */
{
  char *str;
  int len, nhrs, nmin, nsec;

/* -- time interval in code units (dt) -- */

  output->dt = atof(ParamFileGet(output_format, 1));

/* -- check the 2nd field and decide to set "dn" or "dclock" -- */

  str = ParamFileGet(output_format,2);
  len = strlen(str);
  if (str[len-1] == 'h'){
    output->dclock = atof(str);    /* clock interval in hours */
    nhrs = (int)output->dclock;    /* integer part */
    nmin = (int)((output->dclock - nhrs)*100.0); /* remainder in minutes */
    if (nmin >= 60){
      printf ("! OutputFrequency: number of minutes exceeds 60 in %s output\n",
               output_format);
      QUIT_PLUTO(1);
    }
    output->dclock = nhrs*3600.0 + nmin*60;  /* convert to seconds */
    output->dn     = -1;
  }else if (str[len-1] == 'm'){
    output->dclock = atof(str);      /* clock interval in minutes */
    nmin = (int)output->dclock;      /* integer part */
    nsec = (int)((output->dclock - nmin)*100.0); /* remainder in seconds */
    if (nsec >= 60){
      printf ("! OutputFrequency: number of seconds exceeds 60 in %s output\n",
               output_format);
      QUIT_PLUTO(1);
    }
    output->dclock = nmin*60.0 + nsec;
    output->dn     = -1;
  }else if (str[len-1] == 's'){
    output->dclock = atof(str);           /* clock interval in seconds */
    output->dn     = -1;
  }else{
    output->dclock = -1.0;    
    output->dn     = atoi(ParamFileGet(output_format, 2));
  }

}
