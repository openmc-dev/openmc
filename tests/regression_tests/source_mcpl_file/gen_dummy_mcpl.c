#include <mcpl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int batches=10;
const int particles=10000;

double rand01(){
  double r=rand()/((double) RAND_MAX);
  return r;
}


int main(int argc, char **argv){
  char outfilename[128];
  snprintf(outfilename,127,"source.%d.mcpl",batches);

  /*generate an mcpl_header_of sorts*/
  mcpl_outfile_t outfile=mcpl_create_outfile(outfilename);
  mcpl_particle_t *particle, Particle;

  char line[256];
  snprintf(line,255,"gen_dummy_mcpl.c");
  mcpl_hdr_set_srcname(outfile,line);
  mcpl_enable_universal_pdgcode(outfile,2112);/*all particles are neutrons*/
  snprintf(line,255,"Dummy MCPL-file for testing OpenMC mcpl input");
  mcpl_hdr_add_comment(outfile,line);

  /*also add the instrument file and the command line as blobs*/
  FILE *fp;
  if( (fp=fopen("gen_dummy_mcpl.c","rb"))!=NULL){
    unsigned char *buffer;
    int size,status;
    /*find the file size by seeking to end, "tell" the position, and then go back again*/
    fseek(fp, 0L, SEEK_END);
    size = ftell(fp); // get current file pointer
    fseek(fp, 0L, SEEK_SET); // seek back to beginning of file
    if ( size && (buffer=malloc(size))!=NULL){
      if (size!=(fread(buffer,1,size,fp))){
        fprintf(stderr,"Warning: c source generator file not read cleanly\n");
      }
      mcpl_hdr_add_data(outfile, "mcpl_point_source_file_generator", size, buffer);
      free(buffer);
    }
    fclose(fp);
  } else {
    fprintf(stderr,"Warning: could not open c source generator file, hence not embedded.\n");
  }


  /*the main particle loop*/
  particle=&Particle;
  int i;
  for (i=0;i<batches*particles; i++){
    if(i%particles==0){
      printf("batch number %d\n",i/particles);
    }
    /*generate a point source - i.e. x,y,z are 0*/
    particle->position[0]=0;
    particle->position[1]=0;
    particle->position[2]=0;
    
    /*generate a random direction on unit sphere*/
    double nrm=2.0;
    double vx,vy,vz;
    int iter=0;
    do {
      if(iter>100){
        printf("warning: exceeed max random vector attempts: at loop %d -> %g   %d\n",i,nrm,iter);
      }
      vx=rand01()*2.0-1.0;
      vy=rand01()*2.0-1.0;
      vz=rand01()*2.0-1.0;
      nrm=(vx*vx + vy*vy + vz*vz);
      iter++;
    } while (nrm>1.0);
    vx/=sqrt(nrm);
    vy/=sqrt(nrm);
    vz/=sqrt(nrm);
    particle->direction[0]=vx;
    particle->direction[1]=vy;
    particle->direction[2]=vz;

    /*generate the kinetic energy (in MeV) of the particle*/ 
    particle->ekin=1e-3;
 
    /*set time=0*/
    particle->time=0;
    /*set weight*/
    particle->weight=1;

    particle->userflags=0;
    mcpl_add_particle(outfile,particle);
  }
  mcpl_closeandgzip_outfile(outfile);
}
