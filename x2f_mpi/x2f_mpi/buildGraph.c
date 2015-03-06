#include "graphtypes.h"
#include <stdio.h>
#include <math.h>

Graph buildGraph(ngrp,dgr,nds)
int *ngrp;
int *dgr, *nds;
{
	Graph graph;
    int degree, adj_node, elabel, i, j;
    int *p2,*p3;

	p2 = (int *)dgr;
	p3 = (int *)nds;
	graph =  NewGraph(*ngrp);
	/*printf("ngrp is %d\n:",*ngrp);*/
	for (i=1; i<= *ngrp; ++i)
	{
		NLabel(graph,i) = 3;
		Xcoord(graph,i) = 0;
		Ycoord(graph,i) = 0;
		degree = *(p2++);
		/*printf("degree are %d\n",degree);*/
        for (j=1; j<=degree; ++j)
		{
			adj_node = *(p3++);
			elabel = 0;
			/*printf("adj_node and elabel %d %d\n", adj_node,elabel);*/
			{
                if  (i<adj_node)
                    AddEdge (graph,i,adj_node,elabel);   
			}
		}
	}
    return(graph);
}
