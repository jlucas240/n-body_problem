#include <stdio.h> 
#include <stdbool.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>

/*
********************************************************NOTES FOR USE**********************************************************
	- The width(x) and length(y) should be a power of 2 just to make life simple
	- The representative graph starts at a 1;
*/

struct quad_tree* readInData(char* f, struct quad_tree* q);
double get_force(struct quad_tree *tree);

struct node {
	double mass;
	long x; 
	long y; 
	double x_force; //force in x dimension
	double y_force; //force in y dimension
	double x_velocity; //velocity 
	double y_velocity;
	double x_coordinates; //coordinates
	double y_coordinates;
};

struct quad_tree {
	struct quad_tree* c0;
	struct quad_tree* c1;
	struct quad_tree* c2;
	struct quad_tree* c3;
	struct node* value;
	int low_x;
	int low_y;
	int high_x; // half of width
	int high_y; // length
	int size;
	float cent_M; // center of mass at 
	struct node_list *quadtrees_nodes; // list of structs which hold the number of nodess + info
};

struct node_list
{
    struct node * node_info[5000];// hold list of node structs
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
		if(q->value != NULL)
		{

			printf("colishion, i need a bucket dady.... %ld\n", q->value->x);
			}
		q->value = v;
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


/* using for testing  can change how ever you please*/

void main() {
	struct quad_tree* Q = init_quadTree(-536870911, 536870912, -536870911, 536870912);
	readInData("results.txt", Q); 
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
	Gm, q, inclination, longitude of the ascending node       **** this can be changed at any time.******
*/

//************************* can be parallelized 

bool intiate_simulation() { // might need file loccation as a parameter takes a quadtree pointer

	return false;
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
		insert(q, val); // insert in to table
		//exit(0);
	}

	fclose(inFile);
	return NULL;
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

/**************************************
This function takes in a quadtree struct, goes down the quadtree and finds the force, mass, velocity
calculates new position
***************************************/

struct quad_tree * movement (struct quad_tree *tree)
{

	double x_var [tree->quadtrees_nodes->num_nodes]; //get list of x components
    double y_var [tree->quadtrees_nodes->num_nodes]; //get list of y components
	double x_node_force; 
	double y_node_force;
	double node_mass;
	double x_node_velocity;
	double y_node_velocity;
	int aa;
	int bb;
    
	
	//parallelize for loop
	//loop through each quadtree, get mass, get velocity, get coordinates, calculate force
    	for (aa = 0; aa < tree->quadtrees_nodes->num_nodes; aa++)
    		{
        	node_mass = tree->quadtrees_nodes->node_info[aa]->mass; //get node mass
			x_node_force, y_node_force = get_force(tree); //get force components on node

        	x_node_velocity = tree->quadtrees_nodes->node_info[aa]->x_velocity + x_node_force/node_mass; //get new  x velocity
        	y_node_velocity = tree->quadtrees_nodes->node_info[aa]->y_velocity + y_node_force/node_mass; //get new y velocity
        
        	x_var[aa] = tree->quadtrees_nodes->node_info[aa]->x_coordinates + x_node_velocity; //get new x coordinates
        	y_var[aa] = tree->quadtrees_nodes->node_info[aa]->y_coordinates + y_node_velocity; //get new y coordinates
        
        	tree->quadtrees_nodes->node_info[aa]->x_velocity = x_node_velocity; //update x velocity
       		tree->quadtrees_nodes->node_info[aa]->y_velocity = y_node_velocity;	//update y velocity
		}
    
	//
    	for (bb = 0; bb < tree->quadtrees_nodes->num_nodes; bb++)
   	 	{
        	tree->quadtrees_nodes->node_info[bb]->x_coordinates  = x_var[bb]; //update each quadtree with new coordinate
        	tree->quadtrees_nodes->node_info[bb]->y_coordinates = y_var[bb];
   	 	}
    

	return tree;

}

double get_force(struct quad_tree *tree)
{
	double x_force;
	double y_force;

return (x_force, y_force);


}

// go to bottom of tree then come back up calculating center of mass at each point

//************************* can be parallelized 

float dfs(struct quad_tree* Q) {
	
	/*if (Q->size == 1) {
		Q->cent_M = Q->value.
	}*/
}