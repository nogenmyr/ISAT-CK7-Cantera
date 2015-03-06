#include "graphtypes.h"

void trymatch(ngrp,dgr,nds,MAT)
int *ngrp, *dgr, *nds, *MAT; 

{
	int i;
	Graph graph, buildGraph();
	int *MAT1;
	
	/*printf("reach here1\n");*/
	graph = buildGraph(ngrp,dgr,nds);
	/*printf("reach here2\n");*/
	MAT1 = Match(graph);
	
    for (i=1; i<=*ngrp; ++i)
	{
		/* printf("%d %d\n",i,MAT1[i]); */
		MAT[i-1] = MAT1[i];
	}
	/*printf("reach here3\n");*/
}

void trymatch_ (ngrp,dgr,nds,MAT)
int *ngrp, *dgr, *nds, *MAT; {
trymatch (ngrp,dgr,nds,MAT);
}
