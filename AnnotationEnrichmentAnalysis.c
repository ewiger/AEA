#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <cstring>
#include <time.h>
#include <getopt.h>
using namespace std;

#define CR 13
#define LF 10
#define BOLD "\033[1m"
#define NORMAL "\033[0m"


// file info
FILE *fid;
char index_file[500];
char annotation_file[500];
char signature_file[500];
char partition_file[500];
char out_file[500];
char pname[500];


// counters
int v,c,p;
int q1, q2;
int cnt, cnt1, cnt2, randcnt;


// random variables (used in temporary storage capacity)
char temp[999998];
string stemp;
int ctemp[25000];
int gtemp[20000];
int ttemp[20000];

// vectors to store randomizations
int RandTermSignature[20000];
int RandGeneSignature[20000];
int LocRandVec[20000];


// program defaults
int numrandomizations=10000;
int randomcommunities=0;
int randomsignatures=0;
int evaluation_type=0;
int analytic_approximation=0;
int minannotations=0;
int minprogeny=0;


// values calculates during AEA
int NumGenes, NumTerms, NumAnnotations, NumAnnotatedTerms, NumDegreeGenes, NumSignatures, NumCommunities, NumRemaining;
float meandegree, expdegree;
double p_value[16000][600];
int NumShared[16000][600];


// structures used in AEA

typedef struct{
char ID[10];		// GO term ID
int genes[20000];	// genes annotated to GO term
int num;		// number of genes annoted to term
} ANNOTATIONS;
ANNOTATIONS Annotations[16000];

typedef struct{
char name[100];
int index;
int degree;	// number of terms annotated to gene, read in from file
int num;	// number of terms annotated to gene, determined from annotation file (should be equivalent to degree)
int annotations[1000];	// terms to which gene is annotated
} GENES;
GENES Genes[20000];

typedef struct{
char name[1000];
int genes[1000];	// list of genes in signature
int num;	// number of genes in signature
int degree;	// number of annotations made to signature
} SIGNATURES;
SIGNATURES Signatures[600];


typedef struct{
char name[100];
int terms[12000];	// list of terms in community/branch
int num;	// number of terms in community/branch
int numedges;	// number of annotations made to community/branch
} COMMUNITIES;
COMMUNITIES Communities[16000];


// functions
void useage();
double logsum(int start_value, int end_value);
double pval(int N1, int N2, int Overlap, int Total);
int sort(const void *x, const void *y){return (*(int*)x - *(int*)y);}



int main(int argc, char *argv[])
{
	extern char *optarg;
	int cinput;
	int errflg=5;
	strcpy(pname, argv[0]);

	while((cinput = getopt(argc, argv, "i:o:a:p:s:n:t:r:m:M")) != -1)
	switch (cinput)
	{
		case 'i':
			strcpy(index_file, optarg); errflg--;
			break;
		case 'a':
			strcpy(annotation_file, optarg); errflg--;
			break;
		case 's':
			strcpy(signature_file, optarg); errflg--;
			break;
		case 'p':
			strcpy(partition_file, optarg); errflg--;
			break;
		case 'o':
			strcpy(out_file, optarg); errflg--;
			break;
		case 'n':
			numrandomizations=(int) atof(optarg);
			break;
		case 't':
			evaluation_type=atoi(optarg);
			if(evaluation_type==0) analytic_approximation=1;
			break;
		case 'r':
			randomcommunities=atoi(optarg);
			if(randomcommunities==0) randomsignatures=1;
			if(randomcommunities==2) {randomcommunities=1; randomsignatures=1;}
			break;
		case 'm':
			minprogeny=atoi(optarg);
			break;
		case 'M':
			minannotations=atoi(optarg);
			break;
		default:
			errflg++;
	}

	if(errflg)
	{

		useage();
		exit(2);
	}

	if(evaluation_type!=0 && evaluation_type!=1 && evaluation_type!=2 && evaluation_type!=3)
	{
		printf("ERROR: Wrong value specificied for evalution type.  Please indicate 0 for gene randomization or 2 for term randomization.\n");
		exit(1);
	}
	// else if(evaluation_type==0) analytic_approximation=1;
	// else analytic_approximation=0;

	if(analytic_approximation!=0)
	{
		if(analytic_approximation==1) fprintf(stderr, "Analytic Approximation Selected.  Note this will ignore the -n flag.\n");
		else
		{
			printf("ERROR: Only a value of 0 or 1 can be used to specify usage of the analytic approximation.\n");
			exit(1);
		}
	}

	fprintf(stderr, "Reading in gene information\n");

	if((fid=fopen(index_file,"r"))==NULL)
	{
		printf("ERROR: Unable to open file containing index values.  Please check file name and path and try again.\n");
		exit(1);
	}
	cnt=0;
	
	meandegree=0;
	NumAnnotations=0;
	while(fscanf(fid, "%u\t%99[^\t]%u", &cnt1, temp, &cnt2)==3)
	{
		Genes[cnt1].num=0;
		Genes[cnt1].index=cnt1;
		Genes[cnt1].degree=cnt2;
		strcpy(Genes[cnt1].name, temp);
		NumGenes++;
		meandegree=meandegree+cnt2;
		NumAnnotations=NumAnnotations+cnt2;
		if(cnt2>0) NumDegreeGenes++;
	}
	fclose(fid);
	fprintf(stderr, "Number of genes: %u\n", NumGenes);
	meandegree=meandegree/NumDegreeGenes;

	if((fid=fopen(annotation_file,"r"))==NULL)
	{
		printf("ERROR: Unable to open file containing term annotations.  Please check file name and path and try again.");
		exit(1);
	}

	cnt=0;
	NumAnnotatedTerms=0;
	while(fgets(temp, 999998, fid)!=NULL)
	{
		stemp=temp;
		strcpy(Annotations[cnt].ID,(stemp.substr(0,10)).c_str());
		Annotations[cnt].num=0;
		cnt2=10;
		cnt1=stemp.find(",", cnt2);
		while(cnt1>0 && cnt1>cnt2+1)
		{
			Annotations[cnt].genes[Annotations[cnt].num]=atoi((stemp.substr(cnt2+1,cnt1-cnt2-1)).c_str());
			Genes[Annotations[cnt].genes[Annotations[cnt].num]].annotations[Genes[Annotations[cnt].genes[Annotations[cnt].num]].num]=cnt;
			Genes[Annotations[cnt].genes[Annotations[cnt].num]].num++;
			cnt2=cnt1;
			cnt1=stemp.find(",",cnt2+1);
			Annotations[cnt].num++;
		}
		if(Annotations[cnt].num>0) NumAnnotatedTerms++;
		cnt++;
	}
	fclose(fid);
	NumTerms=cnt;
	
	NumAnnotatedTerms=0;
	expdegree=0;
	for(cnt=0; cnt<NumTerms; cnt++)
	{
		if(Annotations[cnt].num>0)
		{
			NumAnnotatedTerms++;
			expdegree+=Annotations[cnt].num;
		}
	}
	expdegree=expdegree/NumAnnotatedTerms;

	if((fid=fopen(signature_file, "r"))==NULL)
	{
		printf("ERROR: Unable to open file containing gene signatures.  Please check file name and path and try again.");
		exit(1);
	}

	if(randomsignatures==1)
	{
		srand(100);
		for(c=0; c<NumGenes; ++c)
		{
			v= rand() % (c+1);
			RandGeneSignature[c]=RandGeneSignature[v];
			RandGeneSignature[v]=c;
		}
	}

	srand(100);
	cnt=0;
	while(fgets(temp, 99998, fid)!=NULL)
	{
		Signatures[cnt].num=0;
		Signatures[cnt].degree=0;
		stemp = temp;
		cnt2=stemp.find("\t");
		strcpy(Signatures[cnt].name, (stemp.substr(0, cnt2-1)).c_str());
		cnt1=stemp.find(",", cnt2);
		while(cnt1>0 && cnt1>cnt2+1)
		{
			for(c=0; c<NumGenes; c++)
			{
				if(strcmp((stemp.substr(cnt2+1,	cnt1-cnt2-1)).c_str(),Genes[c].name)==0)
				{
					if(randomsignatures==1) c=RandGeneSignature[c];
					if(Genes[c].degree>0)
					{
						Signatures[cnt].genes[Signatures[cnt].num]=Genes[c].index;
						Signatures[cnt].degree=Signatures[cnt].degree+Genes[c].degree;
						Signatures[cnt].num++;
						break;
					}
				}
			}

			cnt2=cnt1;
			cnt1=stemp.find(",",cnt2+1);
		}
		fprintf(stderr, "%u,", Signatures[cnt].num);
		cnt++;
	}
	fclose(fid);
	fprintf(stderr, "Number of Signatures: %d\n", cnt);
	NumSignatures=cnt;

	if((fid=fopen(partition_file,"r"))==NULL)
	{
		printf("ERROR: Unable to open file containing term partitions.  Please check file name and path and try again.");
		exit(1);
	}

	if(randomcommunities==1)
	{
		srand(10);
		for(c=0; c<NumTerms; ++c)
		{
			v= rand() % (c+1);
			RandTermSignature[c]=RandTermSignature[v];
			RandTermSignature[v]=c;
		}
	}

	for(cnt=0; cnt<16000; cnt++) Communities[cnt].num=0;
	for(cnt=0; cnt<16000; cnt++) Communities[cnt].numedges=0;

	NumCommunities=0;
	while(fgets(temp, 999998, fid)!=NULL)
	{
		stemp=temp;
		for(c=0; c<NumTerms; c++)
		{
			if(strcmp((stemp.substr(0,10)).c_str(), Annotations[c].ID)==0)
			{
				cnt=c;
				if(randomcommunities==1) cnt=RandTermSignature[c];
				break;
			}
		}

		cnt2=10;
		cnt1=stemp.find(",", cnt2);

		// this only works for GO branches!!!!!!
		c=atoi((stemp.substr(cnt2+1,cnt1-cnt2-1)).c_str());
		strcpy(Communities[c].name, (stemp.substr(0,10)).c_str());

		while(cnt1>0 && cnt1>cnt2+1)
		{
			c=atoi((stemp.substr(cnt2+1,cnt1-cnt2-1)).c_str());
			if(c>NumCommunities) NumCommunities=c;
			Communities[c].terms[Communities[c].num]=cnt;
			Communities[c].num++;
			Communities[c].numedges+=Annotations[cnt].num;
			cnt2=cnt1;
			cnt1=stemp.find(",",cnt2+1);
		}
		cnt++;
	}
	fclose(fid);
	NumCommunities++;
	
	// below makes assumption that if number of communities isnt equal to the number of terms, then these aren't branches
	if(NumCommunities!=NumAnnotatedTerms)
	{
		for(cnt=0; cnt<NumCommunities; cnt++)
		{
			sprintf(temp, "Community%u", cnt);
			strcpy(Communities[cnt].name, temp);
		}
	}

	fprintf(stderr, "Average Gene Degree: %f\n", meandegree);
	fprintf(stderr, "Average Term Degree: %f\n", expdegree);
	fprintf(stderr, "NumGenes: %u\n", NumGenes);
	fprintf(stderr, "NumTerms: %u\n", NumTerms);
	fprintf(stderr, "NumDegreeGenes: %u\n", NumDegreeGenes);
	fprintf(stderr, "NumDegreeTerms: %u\n", NumAnnotatedTerms);
	fprintf(stderr, "NumAnnotations: %u\n", NumAnnotations);
	fprintf(stderr, "NumCommunities:%u\n", NumCommunities);	
	
	if(analytic_approximation==1) fprintf(stderr, "Running Analytic Approximation to Annotation Enrichment Analysis.\n");
	p=0;
	for(cnt=0; cnt<NumCommunities; cnt++)
	{
		if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
		{
			p++;
			for(cnt2=0; cnt2<NumGenes; cnt2++) ctemp[cnt2]=0;
			for(cnt2=0; cnt2<Communities[cnt].num; cnt2++)
			{
				for(cnt1=0; cnt1<Annotations[Communities[cnt].terms[cnt2]].num; cnt1++)
				{
					ctemp[Annotations[Communities[cnt].terms[cnt2]].genes[cnt1]]++;
				}
			}
			for(cnt1=0; cnt1<NumSignatures; cnt1++)
			{
				NumShared[cnt][cnt1]=0;
				for(cnt2=0; cnt2<Signatures[cnt1].num; cnt2++)
				{
					if(evaluation_type==3){if(ctemp[Signatures[cnt1].genes[cnt2]]>0) NumShared[cnt][cnt1]++;}
					else NumShared[cnt][cnt1]+=ctemp[Signatures[cnt1].genes[cnt2]];
				}
				
				if(analytic_approximation==0)
				{
					if(NumShared[cnt][cnt1]>0) p_value[cnt][cnt1]=0;
					else p_value[cnt][cnt1]=numrandomizations;
				}
				else if(analytic_approximation==1)
				{
					if(NumShared[cnt][cnt1]>0) p_value[cnt][cnt1]=pval(Communities[cnt].numedges, Signatures[cnt1].degree, NumShared[cnt][cnt1], NumAnnotations);
					else p_value[cnt][cnt1]=1;
				}
			}
		}
	}
	if(analytic_approximation==0) fprintf(stderr, "Number of Communities Evaluated:%u\nCalculating P-values Based on %u Randomizations!!\n", p, numrandomizations);
	else fprintf(stderr, "Number of Communities Evaluated:%u\n", p);

	time_t seconds;
	seconds=time(NULL);
	srand(seconds);
	if(analytic_approximation==0 && evaluation_type==0)
	{
		fprintf(stderr, "Randomizing!!\n");
		for(randcnt=0; randcnt<numrandomizations; randcnt++)
		{
			fprintf(stderr, "%u,", randcnt);
			for(c=0; c<NumGenes; c++) RandGeneSignature[c]=0;
			for(c=0; c<NumGenes; ++c)
			{
				v= rand() % (c+1);
				RandGeneSignature[c]=RandGeneSignature[v];
				RandGeneSignature[v]=c;
			}
			
			for(c=0; c<NumTerms; c++) RandTermSignature[c]=0;
			for(c=0; c<NumTerms; ++c)
			{
				v= rand() % (c+1);
				RandTermSignature[c]=RandTermSignature[v];
				RandTermSignature[v]=c;
			}

			for(cnt1=0; cnt1<NumSignatures; cnt1++)
			{
				for(cnt2=0; cnt2<NumTerms; cnt2++) ctemp[cnt2]=0;
				for(cnt2=0; cnt2<NumTerms; cnt2++) ttemp[cnt2]=0;
				for(cnt2=0; cnt2<NumGenes; cnt2++) gtemp[cnt2]=0;
				
				v=0;
				for(cnt2=0; cnt2<NumGenes; cnt2++)
				{
					q2=cnt2+1;
					v+=Genes[RandGeneSignature[cnt2]].degree;
					if(v>Signatures[cnt1].degree)
					{
					
						if((v-Signatures[cnt1].degree)>(Signatures[cnt1].degree-(v-Genes[RandGeneSignature[cnt2]].degree)))
						{
							for(c=0; c<Genes[RandGeneSignature[cnt2]].degree; c++)
							{
								ctemp[Genes[RandGeneSignature[cnt2]].annotations[c]]++;
							}
						}						
						break;
					}
					else
					{	
						for(c=0; c<Genes[RandGeneSignature[cnt2]].degree; c++)
						{
							ctemp[Genes[RandGeneSignature[cnt2]].annotations[c]]++;
						}
					}
				}
				
				for(cnt=0; cnt<NumCommunities; cnt++)
				{
					if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
					{
						if(NumShared[cnt][cnt1]>0)
						{
							p=0;
							v=0;
							for(cnt2=0; cnt2<NumTerms; cnt2++)
							{
								v+=Annotations[RandTermSignature[cnt2]].num;
								if(v>Communities[cnt].numedges)
								{
									if((v-Communities[cnt].numedges)>(Communities[cnt].numedges-(v-Annotations[RandTermSignature[cnt2]].num)))
									{
										p+=ctemp[RandTermSignature[cnt2]];
									}
									break;
								}
								else
								{
									p+=ctemp[RandTermSignature[cnt2]];
								}
							}
							if(p>=NumShared[cnt][cnt1]) p_value[cnt][cnt1]+=1.0;
						}
					}
				}
			}
		}
	}
	else if(analytic_approximation==0 && evaluation_type==1)
	{
		fprintf(stderr, "Randomizing Genes!!\n");
		for(randcnt=0; randcnt<numrandomizations; randcnt++)
		{
			fprintf(stderr, "%u,", randcnt);
			for(c=0; c<NumGenes; c++) RandGeneSignature[c]=0;

			for(c=0; c<NumGenes; ++c)
			{
				v= rand() % (c+1);
				RandGeneSignature[c]=RandGeneSignature[v];
				RandGeneSignature[v]=c;
			}
			for(cnt1=0; cnt1<NumSignatures; cnt1++)
			{
				for(cnt2=0; cnt2<NumTerms; cnt2++) ctemp[cnt2]=0;
				v=0;
				for(cnt2=0; cnt2<NumGenes; cnt2++)
				{
					v+=Genes[RandGeneSignature[cnt2]].degree;			
					if(v>Signatures[cnt1].degree)
					{
					
						if((v-Signatures[cnt1].degree)>(Signatures[cnt1].degree-(v-Genes[RandGeneSignature[cnt2]].degree)))
						{
							for(c=0; c<Genes[RandGeneSignature[cnt2]].degree; c++)
							{
								ctemp[Genes[RandGeneSignature[cnt2]].annotations[c]]++;
							}
						}						
						break;
					}
					else
					{	
						for(c=0; c<Genes[RandGeneSignature[cnt2]].degree; c++)
						{
							ctemp[Genes[RandGeneSignature[cnt2]].annotations[c]]++;
						}
					}
					
					
				}
				
				for(cnt=0; cnt<NumCommunities; cnt++)
				{
					if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
					{
						if(NumShared[cnt][cnt1]>0)
						{
							p=0;
							for(cnt2=0; cnt2<Communities[cnt].num; cnt2++)
							{
								p+=ctemp[Communities[cnt].terms[cnt2]];
							}
							if(p>=NumShared[cnt][cnt1]) p_value[cnt][cnt1]+=1.0;
						}
					}
				}
			}
		}
	}
	else if (analytic_approximation==0 && evaluation_type==2)
	{
		fprintf(stderr, "Randomizing Terms!!\n");
		for(randcnt=0; randcnt<numrandomizations; randcnt++)
		{
			fprintf(stderr, "%u,", randcnt);
			
			for(c=0; c<NumTerms; c++) RandTermSignature[c]=0;
			
			for(c=0; c<NumTerms; ++c)
			{
				v= rand() % (c+1);
				RandTermSignature[c]=RandTermSignature[v];
				RandTermSignature[v]=c;
			}

			for(cnt=0; cnt<NumCommunities; cnt++)
			{
				if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
				{
					for(cnt2=0; cnt2<NumGenes; cnt2++) ctemp[cnt2]=0;
					v=0;
					for(cnt2=0; cnt2<NumTerms; cnt2++)
					{
						v+=Annotations[RandTermSignature[cnt2]].num;					
						if(v>Communities[cnt].numedges)
						{
							if((v-Communities[cnt].numedges)>(Communities[cnt].numedges-(v-Annotations[RandTermSignature[cnt2]].num)))
							{
								for(c=0; c<Annotations[RandTermSignature[cnt2]].num; c++)
								{
									ctemp[Annotations[RandTermSignature[cnt2]].genes[c]]++;
								}
							}
							break;
						}
						else
						{
							for(c=0; c<Annotations[RandTermSignature[cnt2]].num; c++)
							{
								ctemp[Annotations[RandTermSignature[cnt2]].genes[c]]++;
							}
						}

						
					}
					
					for(cnt1=0; cnt1<NumSignatures; cnt1++)
					{
						if(NumShared[cnt][cnt1]>0)
						{
							p=0;
							for(c=0; c<Signatures[cnt1].num; c++)
							{
								p+=ctemp[Signatures[cnt1].genes[c]];
							}
							if(p>=NumShared[cnt][cnt1]) p_value[cnt][cnt1]+=1.0;
						}
					}
				}
			}
		}
	}
	else if (analytic_approximation==0 && evaluation_type==3)
	{
		fprintf(stderr, "Randomizing Genes!! (evaluating gene overlap)\n");
		for(cnt=0; cnt<NumCommunities; cnt++)
		{
			// fprintf(stderr, "Community#%u,", cnt);
			if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
			{
				fprintf(stderr, "%s,\n", Communities[cnt].name);
				for(cnt2=0; cnt2<NumGenes; cnt2++) ctemp[cnt2]=0;
				for(cnt2=0; cnt2<Communities[cnt].num; cnt2++)
				{
					for(cnt1=0; cnt1<Annotations[Communities[cnt].terms[cnt2]].num; cnt1++)
					{
						ctemp[Annotations[Communities[cnt].terms[cnt2]].genes[cnt1]]++;
					}
				}

				for(randcnt=0; randcnt<numrandomizations; randcnt++)
				{
					for(c=0; c<NumGenes; c++) RandGeneSignature[c]=0;

					for(c=0; c<NumGenes; ++c)
					{
						v= rand() % (c+1);
						RandGeneSignature[c]=RandGeneSignature[v];
						RandGeneSignature[v]=c;
					}
					for(cnt1=0; cnt1<NumSignatures; cnt1++)
					{
						if(NumShared[cnt][cnt1]>0)
						{
							v=0;
							p=0;
							for(cnt2=0; cnt2<NumGenes; cnt2++)
							{
								v+=Genes[RandGeneSignature[cnt2]].degree;			
								if(v>Signatures[cnt1].degree)
								{
									if((v-Signatures[cnt1].degree)>(Signatures[cnt1].degree-(v-Genes[RandGeneSignature[cnt2]].degree)))
									{
										if(ctemp[RandGeneSignature[cnt2]]>0) p++;
									}					
									break;
								}
								else
								{	
									if(ctemp[RandGeneSignature[cnt2]]>0) p++;
								}
							}
							if(p>=NumShared[cnt][cnt1]) p_value[cnt][cnt1]+=1.0;		
						}
					}
				}
			}
		}

	}
	fprintf(stderr, "\n");
	
	fid=fopen(out_file, "w");
	for(cnt=0; cnt<NumCommunities; cnt++)
		if(Communities[cnt].numedges>=minannotations && Communities[cnt].num>=minprogeny)
			for(cnt1=0; cnt1<NumSignatures; cnt1++)
			{
				if(analytic_approximation==0)
					fprintf(fid,"%u\t%u\t%s\t%s\t%u\t%u\t%u\t%g\n",cnt1,cnt,Signatures[cnt1].name,Communities[cnt].name, Communities[cnt].numedges,Signatures[cnt1].degree,NumShared[cnt][cnt1],p_value[cnt][cnt1]/numrandomizations);
				else
					fprintf(fid,"%u\t%u\t%s\t%s\t%u\t%u\t%u\t%g\n",cnt1,cnt,Signatures[cnt1].name,Communities[cnt].name, Communities[cnt].numedges,Signatures[cnt1].degree,NumShared[cnt][cnt1],p_value[cnt][cnt1]);
			}
	fclose(fid);
}

void useage()
{
	printf("%s\t-i file name for index values of genes%s\n", BOLD, NORMAL);
	printf("%s\t-a file name for annotations to Term%s\n", BOLD, NORMAL);
	printf("%s\t-p file name for partitioning of Terms into Branches/Communities%s\n", BOLD, NORMAL);
	printf("%s\t-s name of file containing signatures%s\n", BOLD, NORMAL);
	printf("%s\t-o name for output file%s\n", BOLD, NORMAL);
	printf("%s\t-n (optional) number of randomizations (default 10,000)%s\n", BOLD, NORMAL);
	printf("%sMORE ADVANCED OPTIONS (generally leave unspecified, see README for more information):%s\n", BOLD NORMAL);
	printf("%s\t-t (optional) Different analysis type%s\n", BOLD, NORMAL);
	printf("%s\t-r (optional) Randomization options%s\n", BOLD, NORMAL);
}

double pval(int N1, int N2, int Overlap, int Total)
{
	double p_value=0;
	int MostExtremeOverlap;

	if(N1>N2) MostExtremeOverlap=N2;
	else MostExtremeOverlap=N1;

	for(int CurrentOverlap=Overlap; CurrentOverlap<=MostExtremeOverlap; CurrentOverlap++)
		p_value+=exp(logsum(N1-CurrentOverlap+1,N1)+logsum(N2-CurrentOverlap+1,N2)+logsum(Total-N2-N1+CurrentOverlap+1,Total-N2)-logsum(1,CurrentOverlap)-logsum(Total-N1+1,Total));
	return p_value;
}

double logsum(int start_value, int end_value)
{
	double result=0;
	for (int counter=start_value; counter <=end_value; counter++)
		result=result+log(counter);
	return result;
}
