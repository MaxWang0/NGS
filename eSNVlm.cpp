#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <iostream>

#include "geneInfo.cpp"
#include "obsRD.cpp"

#define MAXGENE 300000

int main(int argc, char * argv[]) {
	if(argc!=5) {
		printf("Usage: %s inRDfile alpha odf erd h simuord mani simualpha simuh cnv\n", argv[0]);  // first argument program name
		exit(0);
	}

	int nSample= 0;
	int nSample2 = 0;
	 // sscanf(argv[2], "%d", &nSample);  // manually input the sample number

	ObsReadDepth obsRD(MAXGENE);
	ObsReadDepth obsRD2;

	 nSample2 = obsRD2.getsamplenumber(argv[1]) - 1;

	 // printf("Hello, world! %d \n", nSample2);

	obsRD.inputFromFile(argv[1], nSample2);// read the rpk matrix file

	// printf("Hello, world!\n");

	// printf("Total genes: %d\n", obsRD.get_nGene());
	
//	obsRD.display();
//	obsRD.get_oldMedianSample().display();
//	obsRD.get_medianSample().display();
//	obsRD.get_medianGene().display();

	obsRD.cal_alpha();

	// printf("alpha = %.3lf\n", obsRD.get_alpha().get_data(0));

// display the summary of alpha

    
	 
	 FILE * fpo = fopen(argv[2], "w");
	 if(fpo == NULL) {
		printf("Error: cannot open file %s to save\n", argv[3]);
		exit(1);
	 }

	 // obsRD.get_alpha().display();

	 obsRD.get_alpha().summary(fpo);


	 fclose(fpo);


	 obsRD.cal_odf();




//	obsRD.get_gene(0).display();

      FILE * fpp = fopen(argv[3], "w");
	 if(fpp == NULL) {
		printf("Error: cannot open file %s to save\n", argv[3]);
		exit(1);
	 }

	 for(int i=0; i<obsRD.get_nGene(); i++)
	 	obsRD.get_odf(i).summary(fpp);
	
	  fclose(fpp);

	  FILE * fpq = fopen(argv[4], "w");
	  if(fpq == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[4]);
		  exit(1);
	  }


	  for(int i = 0; i<obsRD.get_nGene(); i++)
		  obsRD.get_erd(i).summary(fpq);
	
	  fclose(fpq);
/**

	  obsRD.cal_h();


	  FILE * fpr = fopen(argv[5], "w");
	  if(fpr == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[5]);
		  exit(1);
	  }

	   for(int i = 0; i<obsRD.get_nGene(); i++)
		   obsRD.get_h(i).summary(fpr);
	
	// printf("%d\n", obsRD.get_nGene());

      fclose(fpr);

	  int nmani = obsRD.cal_simu();

	  FILE * fps = fopen(argv[6], "w");
	  FILE * fpm = fopen(argv[7], "w");
	  if(fps == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[6]);
		  exit(1);
	  }
	  if(fpm == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[7]);
		  exit(1);
	  }

	  for(int i = 0; i<1000; i++)
		  obsRD.get_simuord(i).summary(fps);

      obsRD.get_mani().summary(fpm, nmani);
	

	  fclose(fps);
	  fclose(fpm);

	  obsRD.cal_simualpha();
	  
	  FILE * fpsa = fopen(argv[8], "w");
	  if(fpsa == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[8]);
		  exit(1);
	  }

	  for(int i = 0; i<obsRD.get_nGene(); i++)
		  obsRD.get_simualpha().summary(fpsa);


      fclose(fpsa);

	  obsRD.cal_simuh();

	  FILE * fpsh = fopen(argv[9], "w");
	  if(fpsh == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[9]);
		  exit(1);
	  }

	  for(int i = 0; i < obsRD.get_nGene(); i++)
		  obsRD.get_simuh(i).summary(fpsh);

	  fclose(fpsh);


	  int ncnv = obsRD.cal_CNV();

	  FILE * fcnv = fopen(argv[10], "w");
	  if(fcnv == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[10]);
		  exit(1);
	  }

	  obsRD.get_cnv().summary(fcnv, ncnv);

	  

	  fclose(fcnv);

	  int nsen = obsRD.map(nmani, ncnv);

	  FILE * fsen = fopen(argv[11], "w");
	  if(fsen == NULL) {
		  printf("Error: cannot open file %s to save\n", argv[11]);
		  exit(1);
	  }

	  obsRD.get_sen().summary()









	
	//	obsRD.get_odf(i).display();

	printf("%.3lf\n", obsRD.get_odf(0).get_mean_data());
	printf("%.3lf\n", obsRD.get_odf(0).get_median_data());
*/

	return 1;
}

