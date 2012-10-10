#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <set>
#include "khmer.h"
#include "hashbits.h"

#define rand32 ((((unsigned long int)rand()&(unsigned long int)255)<<(unsigned long int)24)+(((unsigned long int)rand()&(unsigned long int)255)<<(unsigned long int)16)+(((unsigned long int)rand()&(unsigned long int)255)<<(unsigned long int)8)+((unsigned long int)rand()&(unsigned long int)255))

using namespace khmer;
using namespace std;

unsigned char *F1,*F2;
double factor;
int K;
int S;
int N;
double P;
int L;
unsigned long int MASK;
int localSize;

class tOracle{
public:
	Hashbits *myOracle;
	void setupOracle(void);
	bool isKmerSet(unsigned long long int kmer);
	void setKmer(unsigned long long int kmer);
	int density(void);
	void mapBack(unsigned char *who);
} *oracle;

void addKmer(unsigned long int kmer,unsigned char *to){
	to[(kmer>>3)]|=1<<(kmer&7);
}

bool isKmerSet(unsigned long int kmer,unsigned char *to){
	if(((to[(kmer>>3)]>>(kmer&7))&1)==1)
		return true;
	return false;
}

void removeKmer(unsigned long int kmer,unsigned char *from){
	from[(kmer>>3)]&=255-(1<<(kmer&7));
}

unsigned long int randomBuddy(unsigned long int si){
	switch(rand()&1){
		case 0: return ((si>>2)+((rand()&3)<<((K-1)*2)))&MASK; break;
		case 1: return ((si<<2)+(rand()&3))&MASK; break;
	}
	cout<<"shouldn't be here"<<endl;
	return 0;
}

void addRandomSequence(unsigned long int si,int length,unsigned char *f1,unsigned char *f2){
	//cout<<si<<endl;
	if(isKmerSet(si,f1)){
		if(length>1)
			addRandomSequence(randomBuddy(si),length,f1,f2);
	}
	else{
		addKmer(si,f1);
		addKmer(si,f2);
		oracle->setKmer(si);
		//we have to add a node also to the oracle here later
		if(length>1)
			addRandomSequence(randomBuddy(si),length-1,f1,f2);
	}
}

int density(unsigned char *f){
	int i=0;
	int j;
	for(j=0;j<N;j++)
		if(isKmerSet(j,f)) i++;
	return i;
}

void sizeOfThisComponent(unsigned long int who,unsigned char *table){
	int i;
	unsigned long int next;
//	cout<<localSize<<endl;
	for(i=0;i<4;i++){
		next=((who<<2)+i)&MASK;
		if(isKmerSet(next, table)){
			removeKmer(next, table);
			localSize++;
			sizeOfThisComponent(next,table);
		}
	}
	for(i=0;i<4;i++){
		next=((who>>2)+((rand()&3)<<((K-1)*2)))&MASK;
		if(isKmerSet(next, table)){
			removeKmer(next, table);
			localSize++;
			sizeOfThisComponent(next,table);
		}
	}
}
int sizeOfThisComponentLinear(unsigned long int who,unsigned char *table){
	set<unsigned long int> currentNodes,newNodes;
	set<unsigned long int>::iterator I;
	unsigned long int next,w;
	int i;
	int counter=1;
	newNodes.insert(who);
	removeKmer(who, table);
	do{
		currentNodes.clear();
		currentNodes.swap(newNodes);
		newNodes.clear();
		for(I=currentNodes.begin();I!=currentNodes.end();I++){
			w=*I;
			for(i=0;i<4;i++){
				next=((w<<2)+i)&MASK;
				if(isKmerSet(next, table)){
					removeKmer(next, table);
					newNodes.insert(next);
					counter++;
				}
			}
			for(i=0;i<4;i++){
				next=((w>>2)+(i<<((K-1)*2)))&MASK;
				if(isKmerSet(next, table)){
					removeKmer(next, table);
					newNodes.insert(next);
					counter++;
				}
			}
		}
	}while(newNodes.size()!=0);
	return counter;
}

int main (int argc, char * const argv[]) {
	int i,j,lcc,lccO,R,r;
	unsigned long long int si;
	int toAdd;
	int D1,D2;
	map<int,int> LCC_Dist;
	map<int,int>::iterator LI;

	factor=(double)atof(argv[6]);
	K=atoi(argv[5]);
	N=pow(4.0,(double)K);
	MASK=N-1;
	P=(double)atof(argv[1]);
	L=atoi(argv[2]);
	toAdd=N*P;
	R=atoi(argv[3]);
	F1=(unsigned char*)malloc((N>>3)*sizeof(unsigned char));
	F2=(unsigned char*)malloc((N>>3)*sizeof(unsigned char));
	srand( (unsigned)time( NULL ) );
	
	
	FILE *f=fopen(argv[4],"w+t");
	char fn[1000];
	strcpy(fn,"dist_");
	strcat(fn,argv[4]);
	FILE *D=fopen(fn,"w+t");
	cout<<"graph size: "<<N<<" K: "<<K<<endl;
	for(r=0;r<R;r++){ //loops
		oracle=new tOracle;
		oracle->setupOracle();

		for(i=0;i<N>>3;i++){
			//set all to 0
			F1[i]=0;
			F2[i]=0;
		}
		
		j=N*P;
		while(j>L){
			do{
				si=rand32&MASK;
			}while(isKmerSet(si,F1));
			addRandomSequence(si,L,F1,F2);
			j-=L;
		}
		do{
			si=rand32&MASK;
		}while(isKmerSet(si,F1));
		addRandomSequence(si,j,F1,F2);
		cout<<"filled"<<endl;
		for(i=0;i<N>>3;i++){
			F2[i]=0;
		}
		oracle->mapBack(F2);
		D1=density(F1);
		D2=density(F2);
//		cout<<density(F1)<<" "<<density(F2)<<" "<<oracle->density()<<endl;
//		for(i=0;i<N;i++)
//			printf("%i	%i	%i	%i\n",i,isKmerSet(i, F1),isKmerSet(i, F2),oracle->isKmerSet(i));
		lcc=0;
		lccO=0;
		for(si=0;si<(unsigned long int)N;si++){
			if(isKmerSet(si, F1)){
				removeKmer(si,F1);
				//localSize=1;
				//sizeOfThisComponent(si,F1);
				localSize=sizeOfThisComponentLinear(si,F1);
				//cout<<localSize<<endl;
//				fprintf(D,"%i\n",localSize);
				if(localSize>lcc)
					lcc=localSize;
			}
		}
        //*
		LCC_Dist.clear();
		for(si=0;si<(unsigned long int)N;si++){
			if(isKmerSet(si,F2)){
				removeKmer(si,F2);
				localSize=sizeOfThisComponentLinear(si,F2);
				//fprintf(D,"%i\n",localSize);
				LCC_Dist[localSize]++;
				if(localSize>lccO)
					lccO=localSize;
			}
		}
         //*/
		/*
		for(LI=LCC_Dist.begin();LI!=LCC_Dist.end();LI++){
			fprintf(D,"%i	%i\n",LI->first,LCC_Dist[LI->first]);
		}
		 */
        for(LI=LCC_Dist.begin();LI!=LCC_Dist.end();LI++){
//		fprintf(D,"%i	%i	%i	%i\n",LCC_Dist[1],LCC_Dist[2],LCC_Dist[3],LCC_Dist[4]);
//		printf("%i	%i	%i	%i\n",LCC_Dist[1],LCC_Dist[2],LCC_Dist[3],LCC_Dist[4]);
            fprintf(D,"%i   %i\n",LI->first,LI->second);
//            printf("%i   %i\n",LI->first,LI->second);
        }
		cout<<(float)lcc/((float)N*(float)P)<<" "<<(float)lccO/(double)D2<<" "<<(double)D1/(double)N<<" "<<(double)D2/(double)N<<endl;
		fprintf(f,"%f	%f	%f	%f\n",(float)lcc/((float)N*(float)P),(float)lccO/(double)D2,(double)D1/(double)N,(double)D2/(double)N);
		delete oracle;
	}
	fclose(f);
	fclose(D);
    return 0;
}

void tOracle::setupOracle(void){
	vector<khmer::HashIntoType> sizes;
	sizes.push_back((int)(pow(4.0,K)*factor));
	myOracle=new Hashbits(K, sizes);
}

bool tOracle::isKmerSet(unsigned long long int kmer){
	if(myOracle->get_count(kmer)==1)
		return true;
	return false;
}

void tOracle::setKmer(unsigned long long int kmer){
	myOracle->count(kmer);
}
int tOracle::density(void){
	int i=0;
	int j;
	for(j=0;j<N;j++)
		if(myOracle->get_count((unsigned long long int)j)) i++;
	return i;	
}

void tOracle::mapBack(unsigned char *who){
	int j;
	for(j=0;j<N;j++)
		if(myOracle->get_count((unsigned long long int)j))
			addKmer((unsigned long int)j, who);	
}

