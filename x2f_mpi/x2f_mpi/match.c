/* N-cubed Non-weighted matching */
/* Implementation of H. Gabow's labelling scheme */
/* For an explanation of the algorithm see JACM 23, pp221-34 */
/* Written by Edward Rothberg 6/85 */



#include "graphtypes.h"

struct adj {
	int edge;
	struct adj *next;
	};

static struct adj **adj_list;
static struct adj *space;

static int U, V;   /* U = no. of nodes, V = no. of edges */
static int *MATE;  /* stores the mate of each vertex.  0 if vertex is exposed */
static int *FIRST; /* first non-outer vertex in the path back to the start of */
			/* the search */
static int *LABEL; /* multi-purpose label */
static int *OUTER; /* queue of outer vertices to be searched */

static int *END;   /* stores endpoints of edges */
static int qcount = 0;


int *Match (graph)
Graph graph;
{	int x,y,i,u,v,qptr=0;
	struct adj *p;

	U = Degree(graph,0);
	V = NumEdges(graph);

	/* set up matching data structure */
	Initialize (graph);

	/* start off with a greedy matching */
	Greedy();

	/* search for an augmenting path from each exposed node */
	for (u=1; u<=U; ++u) {
		if (MATE[u] != 0)
			continue;
		OUTER[++qcount] = u;
		LABEL[u] = 0;
		while (qcount != qptr) {  /* while queue is not empty */
			x = OUTER[++qptr];
			p = adj_list[x];
			while (p != NULL) {

				/* y is adjacent to an outer vertex */
				y = (END[p->edge]==x) ? END[p->edge-1] : END[p->edge];

				if (MATE[y] == 0 && y != u) {
					/* found an augmenting path */
					MATE[y] = x;
					REMATCH(x,y);
					qptr = qcount;
					break;
					}
				else if (LABEL[y] >= 0) {
					/* created a blossom */
					DOLABEL(x,y);
					}
				else {
					/* extended the search path */
					v = MATE[y];
					if (LABEL[v] < 0) {
						LABEL[v] = x;
						FIRST[v] = y;
						OUTER[++qcount] = v;
						}
					}
				p = p->next; 
			}
		}

		for (i=0; i<=U; ++i)
			LABEL[i] = -1;
		qcount = 0;
		qptr = 0;
		}
	FreeUp();
	return(MATE);
}


/* take input graph and set up matching data structure.  Matching uses a */
/* structure where each edge has a number, and the adjacency list entry  */
/* for an edge contains that number */

Initialize (graph)
Graph graph;

{	int i, j;
	int allocsize;
	int currentedge,adj_node;
	Edge edge;
	struct adj *p;

	currentedge = U+2;
	END = (int *) malloc((U+2*V+1)*(sizeof(int)));
	adj_list = (struct adj **) malloc((U+1)*(sizeof(struct adj *)));
	for (i=1; i<=U; i++)
		adj_list[i] = NULL;
	space = (struct adj *) malloc(2*V*sizeof(struct adj));
	p = space;
	for (i = 1; i <= U; ++i) {
		edge = FirstEdge(graph,i);
		for (j = 1; j <= Degree(graph,i); ++j) {
			adj_node = EndPoint(edge);
			if (i < adj_node) {
				END[currentedge-1] = i;
				END[currentedge] = adj_node;
				p->edge = currentedge;
				p->next = adj_list[i];
				adj_list[i]=p++;
				p->edge = currentedge;
				p->next = adj_list[adj_node];
				adj_list[adj_node] = p++;
				currentedge += 2;
				}
			edge = NextEdge(edge);
			}
		}

	allocsize = (U+1)*(sizeof(int));
	FIRST = (int *) malloc(allocsize);
	MATE  = (int *) malloc(allocsize);
	LABEL = (int *) malloc(allocsize);
	OUTER = (int *) malloc(allocsize);

	for (i = 0; i <= U; ++i) {
		LABEL[i] = -1;
		FIRST[i] = 0;
		MATE[i] = 0;
		}
}

FreeUp ()
{
	free(adj_list);
	free(space);
	free(END);
	free(FIRST);
	free(LABEL);
	free(OUTER);
}


/* greedy matching.  If a vertex is unmatched, check all adjacent vertices */
/* to see if any of them are also unmatched.  If so, match them. */

Greedy()
{	struct adj *p;
	int i,adj_node;

	for (i=1; i<=U; i++) {
		if (MATE[i]!=0)
			continue;
		p = adj_list[i];
		while (p!=NULL) {
			adj_node = (END[p->edge]!=i) ? 
					END[p->edge] : END[p->edge-1];
			if (MATE[adj_node]==0) {
				MATE[i] = adj_node;
				MATE[adj_node] = i;
				break;
				}
			p = p->next;
			}
		}
}


/* Augment the matching along the augmenting path defined by LABEL */

REMATCH(v,w)
int v,w;

{	int t,x,y;

	t = MATE[v];
	MATE[v] = w;
	if (MATE[t] != v)
		return;
	else if (LABEL[v] <= U) {
		MATE[t] = LABEL[v];
		REMATCH (LABEL[v], t);
		return;
		}
	else {
		x = END[LABEL[v]];
		y = END[LABEL[v]-1];
		REMATCH (x, y);
		REMATCH (y, x);
		return;
		}
}


/* Make all non-outer vertices in the blossom outer */
LabelSub (v,edge,join)
int v,edge,join;

{
	while (v != join) {
		LABEL[v] = edge;
		FIRST[v] = join;
		OUTER[++qcount] = v;
		v = FIRST[LABEL[MATE[v]]];
		}
}


/* return the number of the edge between x and y */
Findedge (x,y)
int x,y;

{	int edge;
	struct adj *p;

	p = adj_list[x];
	while (p != NULL) {
		edge = p->edge;
		if (END[edge] == y || END[edge-1] == y)
			return(edge);
		p = p->next;
		}
	return(0);
}


/* x and y are adjacent and both are outer.  Create a blossom. */

DOLABEL (x,y)
int x,y;

{	int r,s;
	int edge, flag, join;
	register temp;
    
	edge = Findedge (x,y);
	flag = -edge;
	r = FIRST[x];
	s = FIRST[y];
	if (r == s)
		return;
	LABEL[r] = flag;
	LABEL[s] = flag;

	if (s != 0) {
		temp = r;
		r = s;
		s = temp;
		}
    
	r = FIRST[LABEL[MATE[r]]];

	while (LABEL[r] != flag) {
		LABEL[r] = flag;
		if (s != 0) {
			temp = r;
			r = s;
			s = temp;
			}
		r = FIRST[LABEL[MATE[r]]];
		}

	join = r;

	LabelSub(FIRST[x], edge, join);
	LabelSub(FIRST[y], edge, join);

	for (s = 1; s <= qcount; s++)
		if (LABEL[FIRST[OUTER[s]]] > 0)
			FIRST[OUTER[s]] = join;

}

