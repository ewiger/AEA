CC=g++


all: AEA

AEA: AnnotationEnrichmentAnalysis.c
		$(CC) AnnotationEnrichmentAnalysis.c -o AEA
