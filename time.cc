#include "time.h"
#include <iostream>
using namespace std;

string Time::cas_string(time_t casdoc) {
	struct tm *datumCas;
	datumCas=localtime(&casdoc);
	return asctime(datumCas);
}

time_t Time::taketime() {
	time_t sekundy;
	time(&sekundy);
	return sekundy;
}

double Time::delta(time_t t1a, time_t t2a) {
	return (t2a-t1a);
}

void Time::CPUtime() {
	double a=delta(firsttime, secondtime);
	if(a<60) cout<<"CPU time: "<<a<<" s"<<endl<<endl; 
	else if(a>60 && a<3600) 
		cout<<"CPU time: "<<int(a/60)<<" min "
			<<a-60*int(a/60)<<" s"<<endl<<endl; 
	else if(a>3600) 
		cout<<"CPU time: "<<int(a/3600)<<" h "
			<<int(a/60)-60*int(a/3600)<<" min "
			<<a-60*int(a/60)<<" s"<<endl<<endl; 
}	

string cas_ASCII() {
	struct tm *datumCas;
	time_t sekundy;	
	time(&sekundy);
	datumCas=localtime(&sekundy);
	return asctime(datumCas);
}
