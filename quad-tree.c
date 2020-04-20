#include <stdio.h> 
#include <stdbool.h>
#include <stdlib.h> 
#include <math.h>

/*
********************************************************NOTES FOR USE**********************************************************
	- The width(x) and length(y) should be a power of 2 just to make life simple
	- The representative graph starts at a 1;
*/



struct node {
	int mass;
	float radius;
	int x;
	int y;
	int z;
};

struct quad_tree {
	struct quad_tree* c0;
	struct quad_tree* c1;
	struct quad_tree* c2;
	struct quad_tree* c3;
	struct node* value;
	int init_x;
	int init_y;
	int max_x; // width
	int max_y; // length
	int size;
	float cent_M;
};

struct quad_tree* init_quadTree(int x, int y) {
	struct quad_tree* root = (struct quad_tree*) malloc(sizeof(struct quad_tree));
	root->c0 = NULL;
	root->c1 = NULL;
	root->c2 = NULL;
	root->c3 = NULL;
	root->init_x = 0;
	root->init_y = 0;
	root->max_x = x;
	root->max_y= y;
	root->size = x * x ;
	root->value = NULL;
	root->cent_M = 0;

	return(root);
}

/*
	inserts the node v into the quad_tree q. 
	if the insertion is success full returns true
*/

bool insert(struct quad_tree* q, struct node* v) {  // re work to deal with x max being /2 at all times
	// if the size is one the bottom of the tree has been reached add the point there.
	if (q->size == 1) {
		q->value = v;
		return true;
	}

	//in quadrent one
	if (v->x > q->max_x / 2 && v->y > q->max_y / 2) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c0 == NULL) {
			q->c0 = init_quadTree(q->max_x , q->max_y );
			q->c0->init_x = 1 + q->max_x/2;
			q->c0->init_y = 1 + q->max_y/2;
			q->c0->size = q->size / 4;
		}
		return(insert(q->c0, v));
	}

	//in quadrent two
	if (v->x <= q->max_x / 2 && v->y > q->max_y / 2) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c1 == NULL) {
			q->c1 = init_quadTree(q->max_x/2, q->max_y) ;
			q->c1->init_x = q->init_x;
			q->c1->init_y = 1+ q->max_y/2;
			q->c1->size = q->size / 4;
		}
		return(insert(q->c1, v));
	}

	//in quadrent three
	if (v->x <= q->max_x / 2 && v->y <= q->max_y / 2) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c2 == NULL) {
			q->c2 = init_quadTree(q->max_x/2, q->max_y/2);
			q->c2->init_x = q->init_x;
			q->c2->init_y = q->init_y;
			q->c2->size = q->size / 4;
		}
		return(insert(q->c2, v));
	}

	//in quadrent four
	if (v->x > q->max_x / 2 && v->y <= q->max_y / 2) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c3 == NULL) {
			q->c3 = init_quadTree(q->max_x, q->max_y / 2);
			q->c3->init_x =1 + q->max_x/2;
			q->c3->init_y = q->init_y;
			q->c3->size = q->size / 4;
		}
		return(insert(q->c3, v));
	}

}

/*
	This function retreves a node in the treee at the  int x and int y. 
	where q is the quad tree being searched.
*/

struct node* get(struct quad_tree* q, int x, int y) { // add safty for out of bound requests
	// if the function is given a tree that is empty
	if (q == NULL)
		return NULL;
	// if the size is equal to one the element is located there
	if (q->size == 1)
		return q->value;
	// quadrent 1
	if (x > q->max_x / 2 && y > q->max_y / 2) {
		if (q->c0 == NULL)
			return NULL;
		return get(q->c0, x, y);

	}

	// quadrent 2
	else if (x <= q->max_x / 2 && y > q->max_y / 2) {
		if (q->c1 == NULL)
			return NULL;
		return get(q->c1, x, y);
	}

	// quadrent 3
	else if (x <= q->max_x / 2 && y <= q->max_y / 2) {
		if (q->c2 == NULL)
			return NULL;
		return get(q->c2, x, y);
	}

	// quadrent 4
	else if (x > q->max_x / 2 && y <= q->max_y / 2) {
		if (q->c3 == NULL)
			return NULL;
		return get(q->c3, x, y);
	}
}


/* ussing for testing */

void main() {
	struct quad_tree* test = init_quadTree(4, 4);
	struct node* nodeTest = (struct node*) malloc(sizeof(struct node));
	nodeTest->mass = 123;
	nodeTest->radius = 15;
	nodeTest->x = 1;
	nodeTest->y = 1;
	bool hold = insert(test, nodeTest);
	struct node* nodeTest2 = get(test, 1, 2);
	printf("%d\n", nodeTest2->mass);
}

/* 
	read data in x lines at a time. use longitude of acending and mean anomaly to 
	make a cordinate system to put the data into the quad tree. Earth will be (0,0)
	when entering into the quad tree th z will be ignored but stored inthe node value.
	this will portntialy make the aproximation less accurate, but no other way was seen
	to acomplizh this besides this. This may cause the tree to need a bucket bust will 
	need to be tested first. ALSO multiply the values by 1000 to get integers.
	x = q * cos(lonuitude of asending node)
	y = q * sin(lonuitude of asending node)
	z = q * sin(inclination)
	q = perihelion distance
	format from file
	m, r, q, inclination, longitude of the ascending node       **** this can be changed at any time.******
*/

//************************* can be parallelized 

bool intiate_simulation() { // might need file loccation as a parameter takes a quadtree pointer

	return false;
}

//divde by 1000 to get float of the og

//************************* can be parallelized 

float simulate() {
	
	bool hold = intiate_simulation();// send a quad tree pointer

	if (hold == false) {
		// error
	}

	// step 1 travers quadtree and get center of mass at each parent node
	// will be done with DFS

}

// go to bottom of tree then come back up calculating center of mass at each point

//************************* can be parallelized 

float dfs(struct quad_tree* Q) {
	
	/*if (Q->size == 1) {
		Q->cent_M = Q->value.
	}*/
}