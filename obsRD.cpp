
class ObsReadDepth {
	int nGene, nSample;
	double pi[10];
	double prior[10];
	GeneInfo * gene;
	GeneInfo * simugene;
	GeneInfo * simuerd;
	GeneInfo * simuodf;
	GeneInfo medianGene;
	GeneInfo oldMedianSample, medianSample;
	GeneInfo alpha;
	GeneInfo simualpha;
	GeneInfo cnv;
	GeneInfo mani;
	GeneInfo sen;

	GeneInfo * ERD;

	GeneInfo * odf;

	GeneInfo * h;

	GeneInfo * simuh;

	

    public:
	ObsReadDepth() {gene=NULL; nGene=0; nSample=0;}
	ObsReadDepth(int n) {
		nGene=n;
		gene=new GeneInfo[n];
		simugene=new GeneInfo[n];
		simuerd=new GeneInfo[n];
		simuodf=new GeneInfo[n];
		odf=new GeneInfo[n];
		ERD = new GeneInfo[n];
		h = new GeneInfo[n];
		simuh = new GeneInfo[n];
		const float PI = 3.1415926;

		double prior[10] = {0.000634, 0.00211, 0.996, 0.000538, 0.000668, 0.0000357, 0.00000752, 0.00000139, 0.000000361, 0.0000000437} ;
	}

	~ObsReadDepth() {if(nGene!=0) {delete [] gene; delete [] odf;}}

	int resize(int n) {if(n>0) nGene=n; return n;}

	int get_nGene() {return nGene;}
	int get_nSample() {return nSample;}
	GeneInfo & get_gene(int i) {if(i>=0 && i<nGene) return gene[i];}
	GeneInfo & get_odf(int i) {if(i>=0 && i<nGene) return odf[i];}
	GeneInfo & get_erd(int i) {if(i>=0 && i<nGene) return ERD[i];}
	GeneInfo & get_h(int i) {if(i>=0 && i<nGene) return h[i];}
	GeneInfo & get_simuord(int i) {if(i>=0 && i<nGene) return simugene[i];}
	GeneInfo & get_simuh(int i) {if(i>=0 && i<nGene) return simuh[i];}

	GeneInfo & get_oldMedianSample() {return oldMedianSample;}
	GeneInfo & get_medianGene() {return medianGene;}
	GeneInfo & get_medianSample() {return medianSample;}
	GeneInfo & get_alpha() {return alpha;}
	GeneInfo & get_simualpha() {return simualpha;}
	GeneInfo & get_cnv() {return cnv;}
	GeneInfo & get_mani() {return mani;}
	GeneInfo & get_sen() {return sen;}

	int inputFromFile(char * fn, int ns);

	void display() {
	  	printf("#Total Gene %d samples %d\n", nGene, nSample);
		for(int i=0; i<nGene; i++) gene[i].display();
	}

	int getsamplenumber(char * fn );

	int cal_median_sample();
	int cal_median_gene();
	int cal_alpha();
	double cal_alpha(GeneInfo &, GeneInfo &);
	int cal_odf();
	int cal_h();
	int cal_simu();
	int cal_simuh();
	int cal_simualpha();
	double cal_simualpha(GeneInfo &, GeneInfo &);
	int cal_CNV();
	int map(int m, int n);
};


const float PI = 3.1415926;

int ObsReadDepth::cal_median_sample() {  //calculate the MRD
	char id[10]="MRD";
	medianSample.init(id, nSample);
	for(int i=0; i<nSample; i++) {
	  	GeneInfo tmp(nGene);
		for(int j=0; j<nGene; j++) 
			tmp.set_data(j, gene[j].get_data(i));
		medianSample.set_data(i, tmp.get_median_data());
	}
	return 1;
}

int ObsReadDepth::cal_median_gene() {  //calculate gene level MRD
	char id[10]="gMRD";
	medianGene.init(id, nGene);
	for(int i=0; i<nGene; i++) 
		medianGene.set_data(i, gene[i].get_median_data());
	return 1;
}


// read the rpk matrix file and calculate the MRD in both gene and sample level

int ObsReadDepth::getsamplenumber(char * fn){
	FILE * fp = fopen(fn, "r");
	if(fp==NULL) return -1;


	char line[10000];

	int i = 0;

	while(!feof(fp)){
		fscanf(fp, "%[^\n]\n", line);

		 char * pch;

		 pch = strtok(line, "\t");

		 while(pch != NULL){
			i++;
			pch = strtok(NULL, "\t");
		 };

		break;

	}

	return i; //return sample number	


}


int ObsReadDepth::inputFromFile(char * fn, int ns) {
	nSample=ns;
	FILE * fp=fopen(fn, "r");
	if(fp==NULL) return -1;

	char oldId[]="oldMRD";
	oldMedianSample.init(oldId, ns);
	char line[100];
	int ng=0;
	while(!feof(fp)) {
	  	fscanf(fp, "%s", line); //read gene index as string type to line
		if(strcmp(line, "gi")==0) { // head line
			for(int i=0; i<ns; i++) fscanf(fp, "%s", line);
		} else {
			double f;
			gene[ng].init(line, ns);   // initiate the certain gene with title and number of samples
			for(int i=0; i<ns; i++) {
				fscanf(fp, "%lf", &f);  // save the rd value to f variable
				gene[ng].set_data(i, f); 
			}
			ng++;
		}
	};
	fclose(fp);
	for(int i=0; i<ns; i++) oldMedianSample.set_data(i, gene[ng-2].get_data(i)); // Get the MRD
	nGene=ng-2;

	cal_median_sample(); // calculate the MRD
	cal_median_gene();

	return ng;
}

int ObsReadDepth::cal_alpha() {
	char id[]="Alpha";
	alpha.init(id, nGene);

	for(int i=0; i<nGene; i++) {
		double a=cal_alpha(gene[i], medianSample);
		alpha.set_data(i, a);
	}
	return 1;
}

int ObsReadDepth::cal_simualpha() {
	char id[]="Alpha";
	simualpha.init(id, nGene);

	for(int i=0; i<1000; i++) {
		double a=cal_alpha(simugene[i], medianSample);
		simualpha.set_data(i, a);
	}
	return 1;
}
	

double ObsReadDepth::cal_alpha(GeneInfo & y, GeneInfo & x) {
	// fit y = alpha * x
	double xx=0.0, xy=0.0;
	for(int i=0; i<nSample; i++) {
		double x1=x.get_data(i);
		double y1=y.get_data(i);
		xx+=(x1*x1);
		xy+=(x1*y1);
	}
	if(xx!=0) return xy/xx;
	else {
/*
		x.display();
		y.display();
		printf("alpha= %.3lf\n", xy/xx);
*/
		return 0.0;
	}
}

int ObsReadDepth::cal_odf() {
	GeneInfo stdev_sample(nSample);
	GeneInfo stdev_gene(nGene);

	for(int i=0; i<nGene; i++) {
		odf[i].init(nSample);
		ERD[i].init(nSample);
		double a=alpha.get_data(i);  // get alpha of certain gene
		for(int j=0; j<nSample; j++) {
			double erd=a*medianSample.get_data(j);  // get erd of certain gene
			odf[i].set_data(j, 0.0);
			ERD[i].set_data(j, erd);
			double z=gene[i].get_data(j)-erd;  // get ord - erd
			if(a==0 || z==0 || erd==0) {
				printf("Warning: i=%d alpha= %.4f z= %.4f erd= %.4f ord= %.4f\n",
					i, a, z, erd, gene[i].get_data(j));
			} else {
				z=z/sqrt(erd); // get z value
				odf[i].set_data(j, z); // store z value
			}
		}
		stdev_gene.set_data(i, odf[i].get_stdev()); // get stdev of z in certain gene
	}
	
	double mean_stdev_gene=stdev_gene.get_mean_data(); // get mean of stdev in all genes
// printf("mean_stdev_gene= %.4lf\n", mean_stdev_gene);

	for(int i=0; i<nSample; i++) {
		GeneInfo tmp(nGene);
		for(int j=0; j<nGene; j++) 
			tmp.set_data(j, odf[j].get_data(i));
		stdev_sample.set_data(i, tmp.get_stdev());
	}
// stdev_sample.display();  // get stdev of z in certain sample

	for(int i=0; i<nGene; i++) {
		double cg=stdev_gene.get_data(i);
		for(int j=0; j<nSample; j++) {
			double cs=stdev_sample.get_data(j);
			double o=cs*cg/mean_stdev_gene; // get the over dispersion factor
			odf[i].set_data(j, o);
		}
// printf("alpha= %.4lf\n", alpha.get_data(i));

	// odf[i].display();
if(cg>20) {
	// printf("cg= %.4lf\n", cg);
	// printf("alpha= %.4lf\n", alpha.get_data(i));
	// gene[i].display();
	// medianSample.display();
}

	}

	return 1;
}
		
int ObsReadDepth::cal_h() {

	double prior[10] = {0.000634, 0.00211, 0.996, 0.000538, 0.000668, 0.0000357, 0.00000752, 0.00000139, 0.000000361, 0.0000000437} ;

	for(int i = 0; i<nGene; i++) {
		h[i].init(nSample);
		for(int j = 0; j < nSample; j++){
			double ord = gene[i].get_data(j);
			double erd = ERD[i].get_data(j);
			double ODF = odf[i].get_data(j);
			double delta = ODF * sqrt(ERD[i].get_data(j));
			double preh[10];
			double sum = 0;



			 for(int k = 0 ; k < 10; k++){
				preh[k] =  (1 / (delta * sqrt(k * PI))) * exp(-0.5 * pow((ord - erd * (0.5 * k))/(delta * sqrt(k * 0.5)), 2)) * prior[k];
				if (k == 0){
					preh[k] = (1 / (delta * sqrt(0.2 * PI))) * exp(-0.5 * pow((ord - erd/10)/(delta * sqrt(0.1)), 2)) * prior[k];
				}

				sum += preh[k];
			  }
				
			   double hvalue = preh[1]/sum;

			//h[i].set_data(j, (1 / (delta * sqrt(1 * PI))) * exp(-0.5 * pow((ord - erd * (0.5 * 1))/(delta * sqrt(1 * 0.5)), 2)) * prior[1]);
			//h[i].set_data(j, (1 / (delta * sqrt(1 * PI))) * exp(-0.5 * pow((ord - erd * (0.5 * 1))/(delta * sqrt(1 * 0.5)), 2)) * prior[1]);
			h[i].set_data(j, hvalue);

		}
		// h[i].display();
	}
		return 1;
}

int ObsReadDepth::cal_CNV() {
	char id[]="CNV"; 
	cnv.init(id, nGene);
	int n = 0;
	for(int i = 0; i < nGene; i++) {
		for(int j = 0; j < nSample; j++) {
			double H = simuh[i].get_data(j);
			if( H >= 0.65 ) {
				cnv.set_data(n, H, i, j);
				n++;
			}
		}
	}
	return n;
}
				
int ObsReadDepth::map(int m, int n) {     // m is the manipulation CNV number, n is the real CNV number
     int x = 0;
     for( int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++ ){
			if(cnv.get_gi(i) == mani.get_gi(j)){
				if(cnv.get_si(i) == mani.get_si(j)){
					sen.set_data(x, mani.get_data(j), mani.get_gi(j), mani.get_si(j));
					x++;
				}
			}
		 }
	 }
	 return x;
}



int ObsReadDepth::cal_simuh() {

	double prior[10] = {0.000634, 0.00211, 0.996, 0.000538, 0.000668, 0.0000357, 0.00000752, 0.00000139, 0.000000361, 0.0000000437} ;

	for(int i = 0; i<nGene; i++) {
		simuh[i].init(nSample);
		for(int j = 0; j < nSample; j++){
			double ord = simugene[i].get_data(j);
			double erd = simuerd[i].get_data(j);
			double ODF = simuodf[i].get_data(j);
			double delta = ODF * sqrt(simuerd[i].get_data(j));
			double preh[10];
			double sum = 0;

			for(int k = 0 ; k < 10; k++){
				preh[k] =  (1 / (delta * sqrt(k * PI))) * exp(-0.5 * pow((ord - erd * (0.5 * k))/(delta * sqrt(k * 0.5)), 2)) * prior[k];
				if (k == 0){
					preh[k] = (1 / (delta * sqrt(0.2 * PI))) * exp(-0.5 * pow((ord - erd/10)/(delta * sqrt(0.1)), 2)) * prior[k];
				}

				sum += preh[k];
			}

			double hvalue = preh[1]/sum;

			simuh[i].set_data(j, hvalue);

		}
	}

	return 1;
	
}


int ObsReadDepth::cal_simu() {

	char id[] = "mani_CNV";
	mani.init(id, nGene);

	const gsl_rng_type * R;

	gsl_rng * r;
	gsl_rng_env_setup();
	R = gsl_rng_default;
	r = gsl_rng_alloc(R);
	// gsl_rng_set (r, time(0));
	int n = 0;

	for(int i = 0; i < 1000; i++){

		int geneindex = rand() % nGene;    //sampleindex in the range of 0 to nGene

		simugene[i].init(nSample);
		simuodf[i].init(nSample);
		simuerd[i].init(nSample);

		for(int j = 0; j < nSample; j++){
			
			double erd = ERD[geneindex].get_data(j);
			double ODF = odf[geneindex].get_data(j);
			double dev = sqrt(erd); 
			double ord = erd + gsl_ran_gaussian(r, dev * ODF);
			
			while ( ord < 0 ){
				
				ord = erd + gsl_ran_gaussian(r, dev * ODF);
				
				}

				simugene[i].set_data(j, ord);
				simuodf[i].set_data(j, ODF);
				simuerd[i].set_data(j, erd);
			}
			int samplei[20];
			srand(time(NULL));
			for (int k = 0; k < 20; k++){
				samplei[k] = rand() % nSample;
				double SIMUORD = simugene[i].get_data(samplei[k]);
				simugene[i].set_data(samplei[k], SIMUORD/2);
				mani.set_data(n, SIMUORD/2, i, samplei[k]);
				n++;
			}

			



				
		}
		
		return n;

}
