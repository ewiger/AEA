#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string>
using namespace std;

#define CR 13
#define LF 10
#define BOLD "\033[1m"
#define NORMAL "\033[0m"

FILE *fid;

char ontology_file[500];
char species_file[500];
char network_file[500];
char index_file[500];
char full_file[500];
char signature_file[500];
char branch[1];

// random variables
char temp[9998];
string stemp;
int* vtemp = new int[2000];
int cnt, cnt1, cnt2;

// variables for BFS
int rootnode;
int Qstart;
int Qend;
int v,c;

int NumGenes, NumTerms, NumAnnotatedTerms;

typedef struct{
char ID[10];
int ontology;
char name[500];
char def[5000];
int isa[50];
int ancestors[200];
int degree;
int numparents;
int numancestors;
int obsolete;
int reference;
} ONTOLOGY;

ONTOLOGY Ontology[200000];

int MaxOnt=150000;

typedef struct{
string name;
int llannotations[1000];
int allannotations[5000];
int lldegree;
int alldegree;
} BIPARTITE;

BIPARTITE Bipartite[25000];

typedef struct{
char ID[10];
int genes[20000];
int num;
} ANNOTATIONS;

ANNOTATIONS Annotations[20000];

int sort(const void *x, const void *y){return (*(int*)x - *(int*)y);}

char pname[500];
void useage();

int main(int argc, char *argv[])
{
	extern char *optarg;
	int cinput;
	int errflg=4;
	strcpy(pname, argv[0]);

	while((cinput = getopt(argc, argv, "o:b:s:i:")) != -1)
	switch (cinput)
	{

		case 'o':
			strcpy(ontology_file, optarg); errflg--;
			break;
		case 'b':
			strcpy(branch, optarg); errflg--;
			break;
		case 's':
			strcpy(species_file, optarg); errflg--;
			break;
		case 'i':
			strcpy(index_file, optarg); errflg--;
			break;
		default:
			errflg++;
	}

	if(errflg)
	{

		useage();
		exit(2);
	}

	if(strcmp(branch, "F")!=0 && strcmp(branch, "P")!=0 && strcmp(branch, "C")!=0 && strcmp(branch, "A")!=0)
	{
		printf("ERROR: Wrong specification for ontology branch.  Please use one of the following:\nA (All)\nP (Biological Process)\nF (Molecular Function)\nC (Cellular Component)");
		exit(1);
	}


	if((fid=fopen(ontology_file,"r"))==NULL)
	{
		printf("ERROR: Unable to open ontology file.  Please check file name and path and try again.");
		exit(1);
	}

	fprintf(stderr, "Reading in Ontology\n");
	NumTerms=0;
	cnt=0;
	while(!feof(fid))
	{
		fgets(temp, 9998, fid);
		stemp = temp;

		if(strncmp(temp, "id:", 3)==0) // record GO ID of term and set cnt to numeric part of ID, initialize degree and obsolete
		{
			cnt=atoi((stemp.substr(7,7)).c_str());
			if(cnt>MaxOnt)
			{
				cnt=cnt % 2000000 + MaxOnt;
			}
				
			strcpy(Ontology[cnt].ID,(stemp.substr(4,10)).c_str());
			Ontology[cnt].numparents=0;
			Ontology[cnt].degree=0;
			Ontology[cnt].obsolete=0;
		}
		else if(strncmp(temp, "name:", 5)==0) // record full name of term
		{
			c=stemp.find("\n", 6);
			strcpy(Ontology[cnt].name,(stemp.substr(6,c-1)).c_str());
		}
		else if(strncmp(temp, "def:", 4)==0) // record definition of term
		{
			c=stemp.find("\n", 5);
			strcpy(Ontology[cnt].def,(stemp.substr(5,c-1)).c_str());
		}
		else if(strncmp(temp, "namespace:", 10)==0) // record ontology: BP=0, MF=1, CC=2
		{
			if(strncmp(temp, "namespace: biological_process", 29)==0)
				Ontology[cnt].ontology=0;
			else if(strncmp(temp, "namespace: molecular_function", 29)==0)
					Ontology[cnt].ontology=1;
			else if(strncmp(temp, "namespace: cellular_component", 29)==0)
					Ontology[cnt].ontology=2;
			else fprintf (stderr, "problem");
		}
		else if(strncmp(temp, "is_a:", 5)==0) // propagate up "is a" relationships
		{
			if(atoi((stemp.substr(9,7)).c_str())<=MaxOnt)
			{
				Ontology[cnt].isa[Ontology[cnt].numparents]=atoi((stemp.substr(9,7)).c_str());
			}
			else
			{
				Ontology[cnt].isa[Ontology[cnt].numparents]=atoi((stemp.substr(9,7)).c_str()) % 2000000 + MaxOnt;
			}
			Ontology[cnt].numparents++;
		}
		else if(strncmp(temp, "relationship: part_of", 21)==0) // propagate up "part of" relationships
		{
			if(atoi((stemp.substr(25,7)).c_str())<=MaxOnt)
			{
				Ontology[cnt].isa[Ontology[cnt].numparents]=atoi((stemp.substr(25,7)).c_str());
			}
			else
			{
				Ontology[cnt].isa[Ontology[cnt].numparents]=atoi((stemp.substr(25,7)).c_str()) % 2000000 + MaxOnt;
			}
			Ontology[cnt].numparents++;
		}
		else if(strncmp(temp, "is_obsolete:", 12)==0) // if term is obsolete, record
		{
			Ontology[cnt].obsolete=1;
		}
		if(cnt>NumTerms) NumTerms=cnt;

	}
	fclose(fid);
	fprintf(stderr, "Number of Terms in ontology file: %d\n", NumTerms);

	fprintf(stderr, "Determining Ancestors\n");
	// use "breadth first search" approach to determine all ancestors of each term.  Ancestors are found using relationships stored in Ontology.isa

	for(rootnode=0; rootnode<NumTerms; rootnode++)
	{
		if(Ontology[rootnode].obsolete==0)
		{
			Ontology[rootnode].numancestors=0;
			int Q[NumTerms];
			int visited[NumTerms];
			for(cnt=0; cnt<NumTerms; cnt++)
			{
				Q[cnt]=0;
				visited[cnt]=0;
			}
			Q[0]=rootnode; Qstart=0; Qend=0;
			visited[rootnode]=1;
			while(Qstart<=Qend)
			{
				v=Q[Qstart];
				Qstart=Qstart+1;
				for(cnt=0; cnt<Ontology[v].numparents; cnt++)
				{
					if(visited[Ontology[v].isa[cnt]]==0)
					{
						Qend=Qend+1;
						Q[Qend]=Ontology[v].isa[cnt];
						visited[Ontology[v].isa[cnt]]=1;
						Ontology[rootnode].ancestors[Ontology[rootnode].numancestors]=Ontology[v].isa[cnt];
						Ontology[rootnode].numancestors++;
					}
				}
			}

		}
	}

	if((fid=fopen(species_file,"r"))==NULL)
	{
		printf("ERROR: Unable to open annotation file.  Please check file name and path and try again.");
		exit(1);
	}

	fprintf(stderr, "Reading in gene annotations\n");
	NumGenes=0;
	cnt=-1;
	while(!feof(fid))
	{

		fgets(temp, 9998, fid);

		if(strncmp(temp, "!", 1) !=0)
		{
			stemp = temp;
			cnt1=stemp.find("\t", stemp.find("\t", 0)+1);
			cnt2=stemp.find("\t", cnt1+1);

			if(cnt<0)
			{
				cnt=0;
				NumGenes++;
				Bipartite[cnt].name=stemp.substr(cnt1+1,cnt2-cnt1-1);
				Bipartite[cnt].lldegree=0;
			}
						
			if(strcmp((stemp.substr(cnt1+1, cnt2-cnt1-1)).c_str(), (Bipartite[cnt].name).c_str())!=0)
			{
				// Initialize next gene
				cnt=NumGenes;
				for(c=0; c<NumGenes; c++)
				{
					if(strcmp((stemp.substr(cnt1+1, cnt2-cnt1-1)).c_str(),(Bipartite[c].name).c_str())==0)
					{
						cnt=c;
						break;
					}
				}
				if(cnt==NumGenes)
				{
					NumGenes++;
					Bipartite[cnt].name=stemp.substr(cnt1+1,cnt2-cnt1-1);
					Bipartite[cnt].lldegree=0;
				}
			}

			cnt1=stemp.find("\t", cnt2+1);

			if((cnt1-cnt2)==1)
			{
				if(strcmp(branch,"A")==0)
				{
					if(atoi((stemp.substr(cnt1+4,7)).c_str())<=MaxOnt)
					{
						Bipartite[cnt].llannotations[Bipartite[cnt].lldegree]=atoi((stemp.substr(cnt1+4,7)).c_str());
					}
					else
					{
						Bipartite[cnt].llannotations[Bipartite[cnt].lldegree]=atoi((stemp.substr(cnt1+4,7)).c_str()) % 2000000 + MaxOnt;
					}
					Bipartite[cnt].lldegree++;
				}
				else
				{
					cnt2=stemp.find("\t", cnt1+1);
					cnt2=stemp.find("\t", cnt2+1);
					cnt2=stemp.find("\t", cnt2+1);
					cnt2=stemp.find("\t", cnt2+1);
					if(strcmp(branch,(stemp.substr(cnt2+1, 1)).c_str())==0)
					{
						if(atoi((stemp.substr(cnt1+4,7)).c_str())<=MaxOnt)
						{
							Bipartite[cnt].llannotations[Bipartite[cnt].lldegree]=atoi((stemp.substr(cnt1+4,7)).c_str());
						}
						else
						{
							Bipartite[cnt].llannotations[Bipartite[cnt].lldegree]=atoi((stemp.substr(cnt1+4,7)).c_str()) % 2000000 + MaxOnt;
						}
						Bipartite[cnt].lldegree++;
					}
				}
			}
		}
	}
	fclose(fid);

	fprintf(stderr, "Determining All Annotations\n");
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		if(Bipartite[cnt].lldegree>0)
		{
			// Sort and unique the ll annotations to the current gene

			qsort(Bipartite[cnt].llannotations, Bipartite[cnt].lldegree, sizeof(int),sort);

			c=0;
			for(v=1; v<Bipartite[cnt].lldegree; v++)
			{
				if(Bipartite[cnt].llannotations[v-1] != Bipartite[cnt].llannotations[v])
				{
					Bipartite[cnt].llannotations[c]=Bipartite[cnt].llannotations[v-1];
					c++;
				}
			}
			Bipartite[cnt].llannotations[c]=Bipartite[cnt].llannotations[Bipartite[cnt].lldegree-1];
			Bipartite[cnt].lldegree=c+1;

			// Determine all annotations to the current gene

			Bipartite[cnt].alldegree=0;
			for(v=0; v<Bipartite[cnt].lldegree; v++)
			{
				for(c=0; c<Ontology[Bipartite[cnt].llannotations[v]].numancestors; c++)
				{
					Bipartite[cnt].allannotations[Bipartite[cnt].alldegree+c]=Ontology[Bipartite[cnt].llannotations[v]].ancestors[c];
				}
				Bipartite[cnt].alldegree+=Ontology[Bipartite[cnt].llannotations[v]].numancestors;
			}

			for(v=0; v<Bipartite[cnt].lldegree; v++)
			{
				Bipartite[cnt].allannotations[Bipartite[cnt].alldegree]=Bipartite[cnt].llannotations[v];
				Bipartite[cnt].alldegree++;
			}

			// Sort and unique the annotations to the current gene
			qsort(Bipartite[cnt].allannotations, Bipartite[cnt].alldegree, sizeof(int),sort);

			c=0;
			for(v=1; v<Bipartite[cnt].alldegree; v++)
			{
				if(Bipartite[cnt].allannotations[v-1] != Bipartite[cnt].allannotations[v])
				{
					Bipartite[cnt].allannotations[c]=Bipartite[cnt].allannotations[v-1];
					c++;
				}
									
			}
			Bipartite[cnt].allannotations[c]=Bipartite[cnt].allannotations[Bipartite[cnt].alldegree-1];
			Bipartite[cnt].alldegree=c+1;
		}
	}

	cnt1=0;
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		for(v=0; v<Bipartite[cnt].alldegree; v++)
		{
			if(Ontology[Bipartite[cnt].allannotations[v]].degree==0)
			{
				Ontology[Bipartite[cnt].allannotations[v]].reference=cnt1;
				strcpy(Annotations[cnt1].ID,Ontology[Bipartite[cnt].allannotations[v]].ID);
				cnt1++;
			}
			Annotations[Ontology[Bipartite[cnt].allannotations[v]].reference].genes[Ontology[Bipartite[cnt].allannotations[v]].degree]=cnt;
			Annotations[Ontology[Bipartite[cnt].allannotations[v]].reference].num++;
			Ontology[Bipartite[cnt].allannotations[v]].degree++;
		}
	}
	NumAnnotatedTerms=cnt1;

	fprintf(stderr, "Number of genes annotated in species file: %d\n", NumGenes);
	fprintf(stderr, "Number of terms which have gene annotations: %d\n", NumAnnotatedTerms);
	
	fprintf(stderr, "Printing Complete Gene Annotations\n");
	sprintf(full_file, "%s.annotations", index_file);
	fid=fopen(full_file, "w");
	for(cnt=0; cnt<NumAnnotatedTerms; cnt++)
	{
		fprintf(fid, "%s\t", Annotations[cnt].ID);
		for(cnt1=0; cnt1<Annotations[cnt].num; cnt1++)
		{
			fprintf(fid, "%d,", Annotations[cnt].genes[cnt1]);
		}
		fprintf(fid, "\n");
		
	}
	fclose(fid);
	sprintf(full_file, "%s.index", index_file);
	fid=fopen(full_file, "w");
	for(cnt=0; cnt<NumGenes; cnt++)
		fprintf(fid, "%d\t%s\t%u\n", cnt, (Bipartite[cnt].name).c_str(), Bipartite[cnt].alldegree);
	fclose(fid);
	sprintf(full_file, "%s.description", index_file);
	fid=fopen(full_file, "w");
	for(cnt=0; cnt<NumTerms; cnt++)
	{
		if(Ontology[cnt].obsolete==0 && Ontology[cnt].degree>0)
		{
			// fprintf(fid, "%d\t%s\t%d\t%d\t%s\t%s\n", cnt,Ontology[cnt].ID,Ontology[cnt].degree,Ontology[cnt].ontology, Ontology[cnt].name, Ontology[cnt].def);
			fprintf(fid, "%d\t%s\t%d\t%d\t%s", cnt,Ontology[cnt].ID,Ontology[cnt].degree,Ontology[cnt].ontology, Ontology[cnt].name);
		}
	}
	fclose(fid);
}

void useage()
{
	printf ("%s Useage %s\n", BOLD, pname, NORMAL);
	printf("%s\t-o ontology file (OBO format)%s\n", BOLD, NORMAL);
	printf("%s\t-b branch of ontology used to calculate networks (A,F,P, or C)%s\n", BOLD, NORMAL);
	printf("%s\t-s species annotation file (GAF format)%s\n", BOLD, NORMAL);
	printf("%s\t-i file name for index values of genes%s\n", BOLD, NORMAL);
}

