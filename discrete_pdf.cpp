#include <iostream>
#include <vector>       // vector
#include <stdlib.h> //argc,argv
#include <math.h> //log2f,exp2f
#include <stdio.h>

#define PI			3.14159265358979323846
#define DEG_X		30.86
#define DEG_Y		23.43
#define NUMSTARS	1340

struct discrete_pdf {
	double* pdf;
	int size;
	discrete_pdf() {size=0;pdf=NULL;}
	~discrete_pdf() {free(pdf);}
	void resize(int sz){
		if (sz>size) {
			pdf=(double*)realloc(pdf,sizeof(double)*sz);
			for (;size<sz;size++) pdf[size]=0.0f;
		}
	}
	void add(int idx, double val) {
		if (idx>=size) resize(idx+1);
		pdf[idx]+=val;
	}
	void reassign(int old_idx, int new_idx, double val) {
		add(new_idx,val);
		pdf[old_idx]-=val;
	}
	/* use discrete convolution; fft version is not numerical stable */
	void operator+= (const discrete_pdf& x) {
		int sy=x.size-1;
		int s=size+sy;
		double* newpdf=(double*)malloc(sizeof(double)*s);
		double sum=0.0;
		//convolve x with (*this)
		for(int i=0;i<s;i++) {
			newpdf[i]=0.0f;
			for(int j=(i>sy?i-sy:0);j<size&&j<=i;j++) {
				newpdf[i]+=pdf[j]*x.pdf[i-j];
			}
			sum+=newpdf[i];
		}
		for(int i=0;i<s;i++) newpdf[i]/=sum;//normalize
		for (s--;newpdf[s]==0.0&&s>0;s--);//trim trailing zeros
		size=s+1;
		free(pdf);
		pdf=(double*)realloc(newpdf,sizeof(double)*size);
	}
	void print() {
		for (int i=0;i<size;i++) printf("%f ",pdf[i]);
		printf("\n");
	}
	double cdf(double score) {
		if (score>size) return 1.0;
		if (score<0) return 0.0;
		double val=0.0;
		int i;
		for (i=0;i<score;i++) val+=pdf[i];
		//interpolate the fractional part
		val-=pdf[i-1]*(i-score);
		return val;
	}
};
int main(int argc, char *argv[] ) {
	unsigned char a=255;
	if (a==(unsigned char)(-1))
		printf("%d\n",(unsigned char)(-1));
	discrete_pdf img_pdf;
	img_pdf.add(0,2);
	img_pdf.reassign(0,2,1);
	img_pdf+=img_pdf;
	for(int i=0;i<img_pdf.size;i++) std::cout <<img_pdf.pdf[i]<<" ";
	return 0;
}
