#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define MAX_SIZE 50000


/* Structures used in the Program */

typedef struct{
  double xi;
  double t;
  double hr;
  double hi;
  double htr;
  double hti;
  double er;
  double ei;
} Grid_Element_t;


typedef struct{
  int  indx_[4];
} Grid_NBR_t;




Grid_Element_t  fGRID[MAX_SIZE];
int             fN_GRID    = 0;
int             fInitState = 0;
char            GPD_GridFiles[10][128];
char            GPD_CurrentGridFile[512];
int             GPD_NFiles;

extern "C"{
double  GPD_GetGridElement(int indx,int u);
double  _f_get_linear_interpolation(double x1, double y1, double x2, double y2, double x3);

void     GPD_Load_File(const char *filename,int n_lines);
void     GPD_MakeValueMap();
int      GPD_GetMapRange(int low_hi);
void     GPD_Init(int model = 0);
int      GPD_FindGridFile(const char *file);

double   GPD_Get_GPD(double xi,double t, int gpd);
double   GPD_Get_GPD_Fast(double xi,double t, int gpd);
void     gpd_init_(int *model);
float    gpd_get_cff_(float *f_xi,float *f_t,int *f_gpd);
    Grid_NBR_t  find_xi_neighbours_fast(double xi,double mt);
}
/*============================================*/
/*
   The Routines Are here...
*/


void gpd_init_(int *model){
  int c_model;
  c_model = *model;
  GPD_Init(c_model);
}

float gpd_get_cff_(float *f_xi,float *f_t,int *f_gpd){
  double c_xi  = (double) *f_xi;
  double c_t   = (double) *f_t;
  int    c_gpd = *f_gpd;
  float  tt_gpd   = (float) GPD_Get_GPD_Fast(c_xi,c_t,c_gpd);
  /*  printf("%f %f %d %f \n",c_xi,c_t,c_gpd,tt_gpd); */
  return tt_gpd;
}

void GPD_Init(int model)
{

  int  status = 0;


  sprintf(GPD_GridFiles[0],"_CFF_grid_no_Dterm.data");
  sprintf(GPD_GridFiles[1],"_CFF_grid_with_Dterm.data");

  GPD_NFiles = 2;

  if(model<0||model>GPD_NFiles){
    printf("\n\nERROR !: Error initializing the package with model %d\n",model);
    printf("\tPlease specify model in a range (0,%d)",GPD_NFiles);
  }

  status = GPD_FindGridFile(GPD_GridFiles[model]);
  if(status==0){
    GPD_Load_File(GPD_CurrentGridFile,1000);
    printf("\n=========\nThe GPD Package was successfuly initialized\n");
    fInitState = 1;
  } else {
    printf("\n\nERROR !: Error initializing the package %d\n",model);
    printf("File==>> %s Not Found...\n\n",GPD_GridFiles[model]);
  }
  
}

int     GPD_FindGridFile(const char *file)
{
  char *env_dir;
  printf("GPD_FINDFILE: Looking for file ==>>> %s\n",file);
  printf("\t Checking current directory. \n");
  sprintf(GPD_CurrentGridFile,"%s",file);
  if(access(GPD_CurrentGridFile,R_OK)==0){
      printf("\tFile Found  ==>> %s \n",GPD_CurrentGridFile);
      return 0;
  }else {
    printf("\t File Not Found in Current directory....\n");
  }
  

  printf("Checking in $CLASDVCS_PDF\n");
  env_dir = getenv("CLASDVCS_PDF");
  if(env_dir!=NULL){
    printf("Lenght = %d\n",strlen(env_dir));
    sprintf(GPD_CurrentGridFile,"%s/%s",env_dir,file);
    if(access(GPD_CurrentGridFile,R_OK)==0){
      printf("File Found  ==>> %s \n",GPD_CurrentGridFile);
      return 0;
    }
  } else {
    printf("GPD_FindGridFile: the enviromental variable CLASDVCS_PDF is not defined\n");
  }
  
  
  return -1;
}

void    GPD_MakeValueMap()
{

}

int     GPD_GetMapRange(int low_hi)
{

}



void  GPD_Print_Element(int indx)
{
  printf("index = %d , %f %f H = %f %f ; Ht = %f %f , E = %f %f\n",indx,fGRID[indx].xi,fGRID[indx].t,fGRID[indx].hr,fGRID[indx].hi,fGRID[indx].htr,fGRID[indx].hti,fGRID[indx].er,fGRID[indx].ei);
}

void  GPD_Load_File(const char *filename,int n_lines){
  
  FILE *fp;
  double x,t,hr,hi,htr,hti,er,ei;
  char   header[512];
  char   comment[128];
  int    f_n_lines, a_model;
  double var_p;
  printf("Loading the grid file =>>>> %s\n",filename);

  fp  = fopen(filename,"r");
  
  fgets(header,521,fp);

  printf("HEADER : %s\n",header);
  sscanf(header,"%s %d %d %lf",comment,&f_n_lines,&a_model,&var_p);
  printf("FILE Header Information\n");
  printf("--------------------------------\n");
  printf("Comment =>>> %s\n",comment);
  printf("# lines =>>> %d\n",f_n_lines);
  printf("Model   =>>> %d\n",a_model);

  fN_GRID = f_n_lines;
  for(int i=0;i<f_n_lines;i++){
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&t,&hr,&hi,&htr,&hti,&er,&ei);
    fGRID[i].xi  = x;
    fGRID[i].t   = t;
    fGRID[i].hr  = hr;
    fGRID[i].hi  = hi;
    fGRID[i].htr = htr;
    fGRID[i].hti = hti;
    fGRID[i].er  = er;
    fGRID[i].ei  = ei;
    /*    print_grid_element(i); */
    if(i%1000==0){
      printf(".");
      fflush(stdout);
    }
  }

  printf("\n done.... \n%d lines loded ..\n",f_n_lines);
  fclose(fp);
}

Grid_NBR_t  find_xi_neighbours_fast(double xi,double mt)
{
    Grid_NBR_t  i_r;
    double   div_xi,t_log, div_log_t;
    int      n_xi_,n_mt_, mt_idx,xi_idx;
    n_xi_ = 359;
    n_mt_ = 359;
    div_xi = xi/0.01;
    t_log  = log(mt+1.);
    div_log_t = t_log/0.005;
    mt_idx     = (int) div_log_t;
    xi_idx     = (int) div_xi;
    xi_idx--;
    i_r.indx_[0] =  (xi_idx+1)*359 + mt_idx+1;
    i_r.indx_[1] =  (xi_idx+1)*359 + mt_idx;
    i_r.indx_[2] =  xi_idx*359 + mt_idx;
    i_r.indx_[3] =  xi_idx*359 + mt_idx+1;
    return i_r;
}


Grid_NBR_t  find_xi_neighbours( double xi, double t )
{
  Grid_NBR_t  i_r;
  double max_diff   = 1000;
  double max_diff_t = 1000;
  int    indx_up    = -1;
  int    i;

  for(i=0;i<fN_GRID;i++){
    /*    printf("%f %f\n",fGRID[i].t-t,fGRID[i].xi-xi); */
    if(fGRID[i].xi-xi>0&&fGRID[i].xi-xi<=max_diff&&
       fGRID[i].t-t>0&&fGRID[i].t-t<=max_diff_t) {
      max_diff = fGRID[i].xi-xi;
      max_diff_t = fGRID[i].t-t;
      indx_up  = i;
    }  
  }

  i_r.indx_[0] = indx_up;

  max_diff   = 1000;
  max_diff_t = -1000;

  for(i=0;i<fN_GRID;i++){
    //    printf("%f %f\n",fGRID[i].t-t,fGRID[i].xi-xi);
    if((fGRID[i].xi-xi>0&&fGRID[i].xi-xi<=max_diff)&&
       (fGRID[i].t-t<0&&fGRID[i].t-t>=max_diff_t)) {
      max_diff   = fGRID[i].xi-xi;
      max_diff_t = fGRID[i].t-t;
      indx_up  = i;
    }  
  }

  i_r.indx_[1] = indx_up;

  max_diff   = -1000;
  max_diff_t = -1000;

  for(i=0;i<fN_GRID;i++){
    //    printf("%f %f\n",fGRID[i].t-t,fGRID[i].xi-xi);
    if((fGRID[i].xi-xi<0&&fGRID[i].xi-xi>=max_diff)&&
       (fGRID[i].t-t<0&&fGRID[i].t-t>=max_diff_t)) {
      max_diff   = fGRID[i].xi-xi;
      max_diff_t = fGRID[i].t-t;
      indx_up  = i;
    }
  }

  i_r.indx_[2] = indx_up;

  max_diff   = -1000;
  max_diff_t = 1000;

  for(i=0;i<fN_GRID;i++){
    //    printf("%f %f\n",fGRID[i].t-t,fGRID[i].xi-xi);
    if((fGRID[i].xi-xi<0&&fGRID[i].xi-xi>=max_diff)&&
       (fGRID[i].t-t>0&&fGRID[i].t-t<=max_diff_t)) {
      max_diff   = fGRID[i].xi-xi;
      max_diff_t = fGRID[i].t-t;
      indx_up  = i;
    }  
  }
  
  i_r.indx_[3] = indx_up;

  return i_r;
}




Grid_NBR_t  GPD_Find_Neighbours( double xi, double t )
{
  Grid_NBR_t  i_r;
  double max_diff[4]   = {1000, 1000,-1000,-1000};
  double max_diff_t[4] = {1000,-1000,-1000, 1000};
  int    indx_up       = -1;
  int    i;

  for(i=0;i<fN_GRID;i++){
    //    printf("%f %f\n",fGRID[i].t-t,fGRID[i].xi-xi);
    if(fGRID[i].xi-xi>0&&fGRID[i].xi-xi<=max_diff[0]&&
       fGRID[i].t-t>0&&fGRID[i].t-t<=max_diff_t[0]) {
      max_diff[0] = fGRID[i].xi-xi;
      max_diff_t[0] = fGRID[i].t-t;
      i_r.indx_[0]  = i;
    }  

    if((fGRID[i].xi-xi>0&&fGRID[i].xi-xi<=max_diff[1])&&
       (fGRID[i].t-t<0&&fGRID[i].t-t>=max_diff_t[1])) {
      max_diff[1]   = fGRID[i].xi-xi;
      max_diff_t[1] = fGRID[i].t-t;
      i_r.indx_[1]  = i;
    }
    if((fGRID[i].xi-xi<0&&fGRID[i].xi-xi>=max_diff[2])&&
       (fGRID[i].t-t<0&&fGRID[i].t-t>=max_diff_t[2])) {
      max_diff[2]   = fGRID[i].xi-xi;
      max_diff_t[2] = fGRID[i].t-t;
      i_r.indx_[2]  = i;
    }

    if((fGRID[i].xi-xi<0&&fGRID[i].xi-xi>=max_diff[3])&&
       (fGRID[i].t-t>0&&fGRID[i].t-t<=max_diff_t[3])) {
      max_diff[3]   = fGRID[i].xi-xi;
      max_diff_t[3] = fGRID[i].t-t;
      i_r.indx_[3]  = i;
    }     
  }
  
  return i_r;
}

double GPD_Get_GPD_Fast(double xi,double t, int gpd)
{
  Grid_NBR_t  idx;
  double      xi1,xi2,t1,t2;
  double      g1,g2;
  double      gp1,gp2;
  double      result;

  idx =  find_xi_neighbours_fast(xi,t);
  
  xi1 = fGRID[idx.indx_[2]].xi;
  xi2 = fGRID[idx.indx_[0]].xi;
  t1  = fGRID[idx.indx_[1]].t;
  t2  = fGRID[idx.indx_[0]].t;

  gp1 = GPD_GetGridElement(idx.indx_[1],gpd);
  gp2 = GPD_GetGridElement(idx.indx_[0],gpd);
  
  g1 = _f_get_linear_interpolation(t1,gp1,t2,gp2,t);

  gp1 = GPD_GetGridElement(idx.indx_[2],gpd);
  gp2 = GPD_GetGridElement(idx.indx_[3],gpd);
  
  g2 = _f_get_linear_interpolation(t1,gp1,t2,gp2,t);

  result = _f_get_linear_interpolation(xi1,g1,xi2,g2,xi);

  return result;
  
}

double GPD_Get_GPD(double xi,double t, int gpd)
{
  Grid_NBR_t  idx;
  double      xi1,xi2,t1,t2;
  double      g1,g2;
  double      gp1,gp2;
  double      result;

  idx = GPD_Find_Neighbours(xi,t);
  
  xi1 = fGRID[idx.indx_[2]].xi;
  xi2 = fGRID[idx.indx_[0]].xi;
  t1  = fGRID[idx.indx_[1]].t;
  t2  = fGRID[idx.indx_[0]].t;

  gp1 = GPD_GetGridElement(idx.indx_[1],gpd);
  gp2 = GPD_GetGridElement(idx.indx_[0],gpd);
  
  g1 = _f_get_linear_interpolation(t1,gp1,t2,gp2,t);

  gp1 = GPD_GetGridElement(idx.indx_[2],gpd);
  gp2 = GPD_GetGridElement(idx.indx_[3],gpd);
  
  g2 = _f_get_linear_interpolation(t1,gp1,t2,gp2,t);

  result = _f_get_linear_interpolation(xi1,g1,xi2,g2,xi);

  return result;
  
}

double  _f_get_linear_interpolation(double x1, double y1, double x2, double y2, double x3)
{
  double eq1 = (y2-y1)/(x2-x1);
  double y3  = eq1*(x2-x3) + y1;
  return y3;
}

double  GPD_GetGridElement(int indx,int u)
{
  double result = 0.;
  if(indx<0||indx>=fN_GRID){
    printf("GetGridElement: index %d is out of range (0,%d)\n",indx,fN_GRID);
    return 0.;
  }
  if(u==0) result = fGRID[indx].hr;
  if(u==1) result = fGRID[indx].hi;
  if(u==2) result = fGRID[indx].htr;
  if(u==3) result = fGRID[indx].hti;
  if(u==4) result = fGRID[indx].er;
  if(u==5) result = fGRID[indx].ei;

  return result;
}

