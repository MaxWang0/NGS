class GeneInfo {
	char name[32];
	int nSample;
	double * data;
	int * gi;   //gene index
	int * si;   //sample index
    public:
	GeneInfo() {nSample=0; data=NULL; strcpy(name, "NA");}
	~GeneInfo() {if(data!=NULL) free(data);}

	GeneInfo(char * n) {strcpy(name, n);}
	GeneInfo(int n) {nSample=n; data=(double *)calloc(n, sizeof(double)); gi=(int *)calloc(n, sizeof(int)); si=(int *)calloc(n, sizeof(int));}
	GeneInfo(char * id, int n) {
		strcpy(name, id); 
		nSample=n; 
		data=(double *)calloc(n, sizeof(double));
	    gi=(int *)calloc(n, sizeof(int));
		si=(int *)calloc(n, sizeof(int));
	}

	char * get_name() {return name;}
	int get_nSample() {return nSample;}
	double get_data(int i) {if(i>=0 && i<nSample) return data[i];}
	double get_gi(int i) {if(i>=0 && i<nSample) return gi[i];}
	double get_si(int i) {if(i>=0 && i<nSample) return si[i];}

	int set_data(int i, double f) {if(i>=0 && i<nSample) data[i]=f; return 1;}

	int set_data(int i, double f, int g, int h) {if(i>=0 && i<nSample) data[i]=f; gi[i]=g; si[i]=h;  return 1;}

	int init(int n) {nSample=n; data= (double *)calloc(n, sizeof(double)); return 1; }
	
	//(double *)calloc(n, sizeof(double)); return 1;}
	
	int init_name(char * n) {strcpy(name, n); return 1;}
	int init_size(int n) {nSample=n; data=(double *)calloc(n, sizeof(double)); return 1;}
	int init(char * id, int n) {
                strcpy(name, id);
                nSample=n;
                data=(double *)calloc(n, sizeof(double));
				gi=(int *)calloc(n, sizeof(int));
				si=(int *)calloc(n, sizeof(int));
		return 1;
	}

	void display() {
		printf("%-31s", name);
		for(int i=0; i<nSample; i++)
			printf("\t%.1f", data[i]);
		printf("\n");
	}


	double get_median_data();
	double get_mean_data();
	double get_stdev();
	void summary();
	void summary(FILE * );
	void summary(FILE *, int);
};

void GeneInfo::summary(){
	printf("Hello, world!");
}

void GeneInfo::summary(FILE * fp){
	fprintf(fp, "%-31s", name);
	for(int i=0; i<nSample; i++)
		fprintf(fp, "\t%.10f", data[i]);
		fprintf(fp, "\n");
}

void GeneInfo::summary(FILE * fp, int n) {
	fprintf(fp, "%-31s", name);
	for(int i=0; i<n; i++) {
		fprintf(fp, "\t%.3f", data[i]);
		fprintf(fp, "\t%d", gi[i]);
		fprintf(fp, "\t%d", si[i]);
		fprintf(fp, "\n");
	}
}
	

int sort_double(const void * d1, const void * d2) {
 	double dd1=*(const double *)d1;
 	double dd2=*(const double *)d2;
	if(dd1<dd2) return -1;
	else if(dd1>dd2) return 1;
	else return 0;
}

double GeneInfo::get_median_data() {
	if(nSample==0) return 0;
	double * tmp=(double *)calloc(nSample, sizeof(double));
	for(int i=0; i<nSample; i++) tmp[i]=data[i];
	qsort(tmp, nSample, sizeof(double), sort_double);
	if(nSample%2==0) return (tmp[nSample/2-1]+tmp[nSample/2])/2.0;	
	else return tmp[(nSample-1)/2];
}

double GeneInfo::get_mean_data() {
	double sum=0.0;
	for(int i=0; i<nSample; i++) sum+=data[i];
	return sum/nSample;
}

double GeneInfo::get_stdev() {
	if(nSample==0) return 0;
	double sum_x=0.0, sum_x2=0.0;
	for(int i=0; i<nSample; i++) {
		sum_x+=data[i];
		sum_x2+=(data[i]*data[i]);
	}
	double ave=sum_x/nSample;
	double f=sum_x2/nSample-ave*ave;
	return sqrt(f);
}

