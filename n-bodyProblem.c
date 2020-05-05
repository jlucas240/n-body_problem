#include <stdio.h> 
#include <stdbool.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <omp.h>
#include "timer.h"




/*
********************************************************NOTES FOR USE**********************************************************
	- run with gcc  -o n n-bodyProblem.c -lm -g -fopenmp
	- ./n
	
	- The width(x) and length(y) should be a power of 2 just to make life simple
	- The representative graph starts at a 1;
*/

struct quad_tree* readInData(char* f, struct quad_tree* q);
void movement(struct quad_tree* q);
struct node_list *nodesInTree;

struct pair {
	double x;
	double y;
};

typedef struct node {
	double mass;
	long x; 
	long y; 
	double x_force; //force in x dimension
	double y_force; //force in y dimension
	double x_velocity; //velocity 
	double y_velocity;
	double x_coordinates; //coordinates
	double y_coordinates;
} node;

struct pair* get_force(struct quad_tree *q, struct node *n);
struct quad_tree {
	struct quad_tree* c0;
	struct quad_tree* c1;
	struct quad_tree* c2;
	struct quad_tree* c3;
	struct node* value;
	long low_x;
	long low_y;
	long high_x; // half of width
	long high_y; // length
	long size;
	double mass_x;
	double mass_y;
	long cent_M_x; // center of mass at 
	long cent_M_y;
	double mass;
	struct node_list *quadtrees_nodes; // list of structs which hold the number of nodess + info
};

struct node_list
{
    struct node *node_info[960000];// hold list of node structs
    int num_nodes; //number of nodes in list
}; 

// width and length should be the same thing



struct quad_tree* init_quadTree(int low_x, int high_x, int low_y, int high_y) {
	struct quad_tree* root = (struct quad_tree*) malloc(sizeof(struct quad_tree));
	root->c0 = NULL;
	root->c1 = NULL;
	root->c2 = NULL;
	root->c3 = NULL;
	root->low_x = low_x;
	root->low_y = low_y;
	root->high_x = high_x;
	root->high_y = high_y;
	root->size = pow((high_x - low_x+1),2);
	root->value = NULL;
	root->mass_x = 0;	
	root->mass_y = 0;
	root->cent_M_x = 0;
	root->cent_M_y = 0;

	return(root);
}

/*
	inserts the node v into the quad_tree q. 
	if the insertion is success full returns true
*/
omp_lock_t writelock;
//pthread_spinlock_t splock;

bool insert(struct quad_tree* q, struct node* v) {  // re work to deal with x max being /2 at all times
	// if the size is one the bottom of the tree has been reached add the point there.

	if (q->size == 1) {
		if(q->value != NULL)
		{
			printf("colishion, i need a bucket dady.... %ld\n", q->value->x);
			return false;
			}
		//pthread_spin_lock(&splock);
		omp_set_lock(&writelock);
		q->value = v;
		omp_unset_lock(&writelock);
		// pthread_spin_unlock(&splock);
		return true;
	}
	//printf("hold m dick\n");
	double w = sqrt(q->size) / 2;
	//printf("%f\n", w);

	//in quadrent one
	if (v->x > q->high_x-w  && v->y >  q->high_y-w) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c0 == NULL) {
			q->c0 = init_quadTree(q->high_x+1 -w, q->high_x , q->high_y+1-w, q->high_y );
		}
		return(insert(q->c0, v));
	}

	//in quadrent two
	else if (v->x <=  q->high_x-w && v->y >  q->high_y-w) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c1 == NULL) {
			q->c1 = init_quadTree( q->low_x, q->high_x-w,  q->high_y+1 -w , q->high_y);
		}
		return(insert(q->c1, v));
	}

	//in quadrent three
	else if (v->x <= q->high_x - w && v->y <= q->high_y - w) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c2 == NULL) {
			q->c2 = init_quadTree(q->low_x, q->high_x - w, q->low_y, q->high_y - w);
		}
		return(insert(q->c2, v));
	}

	//in quadrent four
	else if (v->x >  q->high_x - w && v->y <=  q->high_y - w) { //does not hava safty for out of bounds
		// if the node is not initiated do so here
		if (q->c3 == NULL) {
			q->c3 = init_quadTree(q->high_x + 1-w, q->high_x, q->low_y, q->high_y-w);
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

	double w = sqrt(q->size) / 2;

	// quadrent 1
	if (x >  q->high_x-w && y > q->high_y-w) {
		if (q->c0 == NULL)
			return NULL;
		return get(q->c0, x, y);

	}

	// quadrent 2
	else if (x <= q->high_x-w && y > q->high_y-w) {
		if (q->c1 == NULL)
			return NULL;
		return get(q->c1, x, y);
	}

	// quadrent 3
	else if (x <=  q->high_x-w && y <=  q->high_y-w) {
		if (q->c2 == NULL)
			return NULL;
		return get(q->c2, x, y);
	}

	// quadrent 4
	else if (x > q->high_x-w && y <= q->high_y-w) {
		if (q->c3 == NULL)
			return NULL;
		return get(q->c3, x, y);
	}
}


/*****************************
 * calculates center of mass and mass at each parent and and leaf
 ******************************/
bool set_up(struct quad_tree *q){

	double c0_cmx = 0, c0_cmy = 0, c0_m = 0;
	double c1_cmx = 0, c1_cmy, c1_m = 0;
	double c2_cmx = 0, c2_cmy, c2_m = 0;
	double c3_cmx = 0, c3_cmy = 0, c3_m = 0;

	if(q->size != 1){
		if(q->c0 != NULL)
			set_up(q->c0);
		if(q->c1 != NULL)
			set_up(q->c1);
		if(q->c2 != NULL)
			set_up(q->c2);
		if(q->c3 != NULL)
			set_up(q->c3);

		if(q->value == NULL){
			if(q->c0 == NULL && q->c1 == NULL && q->c2 == NULL && q->c3 == NULL ){
				return false;
			};
			if(q->c0 != NULL){
				c0_cmx = q->c0->cent_M_x;
				c0_cmy = q->c0->cent_M_y;
				c0_m = q->c0->mass;
			}
			if(q->c1 != NULL){
				c1_cmx = q->c1->cent_M_x;
				c1_cmy = q->c1->cent_M_y;
				c1_m = q->c1->mass;
			}
			if(q->c2 != NULL){
				c2_cmx = q->c2->cent_M_x;
				c2_cmy = q->c2->cent_M_y;
				c2_m = q->c2->mass;
			}
			if(q->c3 != NULL){
				c3_cmx = q->c3->cent_M_x;
				c3_cmy = q->c3->cent_M_y;
				c3_m = q->c3->mass;
			}

			q->mass = c0_m + c1_m + c2_m + c3_m;
			q->cent_M_x = (c0_cmx*c0_m + c1_cmx*c1_m + c2_cmx*c2_m + c3_cmx*c3_m)/q->mass;
			q->cent_M_y = (c0_cmy*c0_m + c1_cmy*c1_m + c2_cmy*c2_m + c3_cmy*c3_m)/q->mass;
			//printf("parent\n");
			return true;
		}
	}

	q->cent_M_x = q->value->x_coordinates;
	q->cent_M_y = q->value->y_coordinates;
	q->mass = q->value->mass;
	//nodesInTree->node_info = q->value;
	return true;

}

bool clear_tree(struct quad_tree *q){


	if(q->c0 != NULL){
		clear_tree(q->c0);
	}
	if(q->c1 != NULL){
		clear_tree(q->c1);
	}
	if(q->c2 != NULL){
		clear_tree(q->c2);
	}
	if(q->c3 != NULL){
		clear_tree(q->c3);
	}

	if(q->value != NULL){
		free(q->value);
	}
	free(q);
	return true;
}


/* using for testing  can change how ever you please*/

void main() {
	struct quad_tree* Q = init_quadTree(-536870911, 536870912, -536870911, 536870912);
	nodesInTree = (struct node_list*) malloc(sizeof(struct node_list));
	double start, finish, elapsed;
	
	//Q = readInData("results.txt", Q); 
	Q = readInData("test.txt", Q);
	set_up(Q);
	GET_TIME(start);
	for (int i = 0; i < 5; ++i){
		movement(Q);
		printf("it worked?\n");
	}
	GET_TIME(finish);
    elapsed = finish - start;
    printf("The code to be timed took %e seconds\n", elapsed);
	
}


// read in data for File f and insert data into quadtree q.
// data formate is 	absolute magnitude parameter, perihelion distance, longitude of the ascending node, semi-major axis

struct quad_tree* readInData(char* f, struct quad_tree* q) {

	char* value = malloc(200);
	FILE* inFile = fopen(f, "r");
	FILE* outx = fopen("outx", "w");
	FILE* outy = fopen("outy", "w");
	int k = 0;

	if (!inFile)
		return NULL;

	while (fgets(value, 150, inFile)) {
		
		int i = 0;
		char * token = strtok(value, ",");
		struct node* val = (struct node*) malloc(sizeof(struct node));
		struct node* val2 = (struct node*) malloc(sizeof(struct node));
		float mass;
		float hold_abs;
		double hold_dis;
		double hold_angle;
		double velocity;
		double a;
		double e;
		double tp;
		double period;
		double r;

		while( token != NULL ) {	
			if(i == 0){ // for absolute magnitude
				hold_abs = strtod(token, NULL);
			}
			else if(i == 1){ // for perihelion distance
				hold_dis = strtod(token, NULL);
			}
			else if(i == 2){ // for longitude of the assending node
				hold_angle = strtod(token, NULL);
			}
			else if(i == 3){ // for semi-major axis
				a = strtod(token, NULL);
			}
			token = strtok(NULL, ",");
			++i;
   		}
		// calculate mass of body
		val->mass = 4.83*(pow(10,(hold_abs-4.83)/-2.5));
		// calculate x position of body
		val->x = (hold_dis*1000000)*cos(hold_angle);
		val->y = (hold_dis*1000000)*sin(hold_angle);
		val->x_coordinates = (hold_dis)*cos(hold_angle);
		val->y_coordinates = (hold_dis)*sin(hold_angle);
		//calculate velocity of body relative to the sun at perihelion distance
		velocity = sqrt(6.67408*pow(10, -11)*1.989*pow(10,30)*((2/(hold_dis*1.496*pow(10,11))-(1/(a*1.496*pow(10,11))))));
		//calculate (x,y) components of velocity
		val->x_velocity = velocity*cos(hold_angle);
		val->y_velocity	= velocity*sin(hold_angle);
		k++;

		if(val->mass == 0)
			printf("WHY");
		
		if(insert(q, val) != false){
			// add node to list of nodes
			nodesInTree->num_nodes = k;
			nodesInTree->node_info[k-1] = val;
		}
		
		 // insert in to tree
	}
	fclose(inFile); 
	return q;
}


/**************************************
This function takes in a quadtree struct, goes down the quadtree and finds the force, mass, velocity
calculates new position
	-this is where parallelism happens
***************************************/

void movement(struct quad_tree* q){
	
	struct node_list *hold_list = (struct node_list*) malloc(sizeof(struct node_list));
	struct quad_tree *hold_tree = init_quadTree(-536870911, 536870912, -536870911, 536870912);
	int j = 0;

	omp_init_lock(&writelock);// lock to provent over writhing data

	// calculations are parallelised
	#pragma omp parallel for num_threads(4) // this valuse can be modified inorder to increase treads**************************
	for (int i  = 0; i < nodesInTree->num_nodes; ++i){
		printf("%d\n",i);
		struct node *hold = nodesInTree->node_info[i]; 
		struct node *hold1 = (struct node*) malloc(sizeof(struct node));

		//if(hold->mass != 0){ // error check
			hold1->mass = hold->mass; // set ne nodes ma
			double x_node_force;
			double y_node_force;
			struct pair* forces = get_force(q, hold);

			// compute force on node
			forces = get_force(q, hold);
			x_node_force = forces->x;
			y_node_force = forces->y;
			free(forces);
			

			//get changes in velocity		
			hold1->x_velocity = x_node_force/hold1->mass;
			hold1->y_velocity =  y_node_force/hold1->mass;

			//get changes in position
			hold1->x_coordinates = hold->x_coordinates+hold1->x_velocity;
			hold1->y_coordinates = hold->y_coordinates+hold1->y_velocity;
			hold1->x = hold1->x_coordinates*1000000;
			hold1->y = hold1->y_coordinates*1000000;

			
			// add to new tree
			if (insert(hold_tree, hold1) != false){
				//add to new list
				hold_list->node_info[j] = hold1;
				hold_list->num_nodes = j + 1;
				++j;
			}				
		//}
		
	}
	
	//relase memory
	for(int i = 0; i < nodesInTree->num_nodes; ++i)
		free(nodesInTree->node_info[i]);
	free(nodesInTree);

	free(q);
	nodesInTree = hold_list;
	q = hold_tree;
	set_up(q);
	omp_destroy_lock(&writelock);
	return;
}

struct pair* get_force(struct quad_tree *q, struct node *n){
	double force_x = 0;
	double force_y = 0;
	struct pair *ret = (struct pair*) malloc(sizeof(struct pair));

	// this when the quad tree block is close enought to the partical or there is no more children
	if(q->size == 1 || q->size/sqrt(pow(q->cent_M_x - n->x_coordinates, 2) + pow(q->cent_M_y - n->y_coordinates, 2)) <= 1){
		force_x = 0.0000000000667*q->mass*(n->x_coordinates- q->cent_M_x)/pow(n->x_coordinates- q->cent_M_x, 3);
		force_y = 0.0000000000667*q->mass*(n->y_coordinates- q->cent_M_y)/pow(n->y_coordinates- q->cent_M_y, 3);
		ret->x = force_x;
		ret->y = force_y;
		return ret;
	}
	struct pair* hold;
	
	
	if(q->c0 != NULL){ // get forces at quadrent 0
		hold = get_force(q->c0,n);
		force_x += hold->x;
		force_y += hold->y;
		free(hold);
	} //
	if(q->c1 != NULL){// get forces at quadrent 1
		hold = get_force(q->c1,n);
		force_x += hold->x;
		force_y += hold->y;
		free(hold);
	}
	if(q->c2 != NULL){// get forces at quadrent 2
		hold = get_force(q->c2,n);
		force_x += hold->x;
		force_y += hold->y;
		free(hold);
	}
	if(q->c3 != NULL){// get forces at quadrent 3
		hold = get_force(q->c3,n);
		force_x += hold->x;
		force_y += hold->y;
		free(hold);
	}
	
	//return the forces
	ret->x = force_x;
	ret->y = force_y;
	return ret;
}
