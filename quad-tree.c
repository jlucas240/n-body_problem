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
void simulat(struct quad_tree* q);
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

	//printf("leaf\n");
	q->cent_M_x = q->value->x_coordinates;
	q->cent_M_y = q->value->y_coordinates;
	q->mass = q->value->mass;
	//nodesInTree->node_info = q->value;
	return true;

}




/* using for testing  can change how ever you please*/

void main() {
	struct quad_tree* Q = init_quadTree(-536870911, 536870912, -536870911, 536870912);
	nodesInTree = (struct node_list*) malloc(sizeof(struct node_list));
	

	
	/*Q = readInData("results.txt", Q); 
	set_up(Q);
	simulat(Q);
	printf("it worked?");*/
	
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
		// add node to list of nodes
		nodesInTree->num_nodes = k;
		nodesInTree->node_info[k-1] = val;
		
		insert(q, val); // insert in to tree
	}
	fclose(inFile); 
	return q;
}


/**************************************
This function takes in a quadtree struct, goes down the quadtree and finds the force, mass, velocity
calculates new position
***************************************/

struct quad_tree * movement (struct quad_tree *tree)
{

	double x_var [nodesInTree->num_nodes]; //get list of x components
    double y_var [nodesInTree->num_nodes]; //get list of y components
	double x_node_force; 
	double y_node_force;
	double node_mass;
	double x_node_velocity;
	double y_node_velocity;
	int aa;
	int bb;
	struct pair* force;
    
	
	//parallelize for loop
	//loop through each quadtree, get mass, get velocity, get coordinates, calculate force
    	for (aa = 0; aa < nodesInTree->num_nodes; aa++)
    		{
        	node_mass = nodesInTree->node_info[aa]->mass; //get node mass
			force = get_force(tree, tree->quadtrees_nodes->node_info[aa]); //get force components on node
			x_node_force = force->x, y_node_force = force->y;

        	x_node_velocity = nodesInTree->node_info[aa]->x_velocity + x_node_force/node_mass; //get new  x velocity
        	y_node_velocity = nodesInTree->node_info[aa]->y_velocity + y_node_force/node_mass; //get new y velocity
        
        	x_var[aa] = nodesInTree->node_info[aa]->x_coordinates + x_node_velocity; //get new x coordinates
        	y_var[aa] = nodesInTree->node_info[aa]->y_coordinates + y_node_velocity; //get new y coordinates
        
        	nodesInTree->node_info[aa]->x_velocity = x_node_velocity; //update x velocity
       		nodesInTree->node_info[aa]->y_velocity = y_node_velocity;	//update y velocity
		}
    
	//
    	for (bb = 0; bb < nodesInTree->num_nodes; bb++)
   	 	{
        	nodesInTree->node_info[bb]->x_coordinates  = x_var[bb]; //update each quadtree with new coordinate
        	nodesInTree->node_info[bb]->y_coordinates = y_var[bb];
   	 	}
    

	return tree;

}

struct pair *get_force(struct quad_tree *tree, struct node *node_n)
{

    double node_mass = node_n->mass;
	double total_mass;
	double length_between_mass;
	double mass_of_x;
	double mass_of_y;


    double node_coordinates_x = node_n->x_coordinates;
    double node_coordinates_y = node_n->y_coordinates;
    double calc_force_x = 0;
    double calc_force_y = 0;
	double x_force_array[4] = {0};
	double y_force_array[4] = {0};

	struct pair *ret = (struct pair*) malloc(sizeof(struct pair)); // structurre used to return 2 elements
    
    if (nodesInTree->num_nodes == 1)
    {
        length_between_mass = sqrt(pow((node_n->x_coordinates - tree->quadtrees_nodes->node_info[0]->x_coordinates), 2) + pow((node_coordinates_y - tree->quadtrees_nodes->node_info[0]->y_coordinates), 2));
        if (length_between_mass > 0)
        {
            calc_force_x = node_mass * tree->mass * 0.00000000006673 * (tree->quadtrees_nodes->node_info[0]->x_coordinates - node_coordinates_x) / pow(length_between_mass, 3);
            calc_force_y = node_mass * tree->mass * 0.00000000006673 * (tree->quadtrees_nodes->node_info[0]->y_coordinates - node_coordinates_y) / pow(length_between_mass, 3);
        }
    }
    else 
    {
       	total_mass = tree->mass;
        mass_of_x = tree->mass_x;
        mass_of_y = tree->mass_y;
        length_between_mass = sqrt(pow((mass_of_x - node_coordinates_x), 2) + pow((mass_of_y - node_coordinates_x), 2));

        if(((tree->high_x - tree->low_x)/length_between_mass) < 0.3 && (length_between_mass > 0) )
        {
            calc_force_x = calc_force_x + (node_mass * total_mass * 0.00000000006673 * (mass_of_x - node_coordinates_x) / pow(length_between_mass, 3));
            calc_force_y =  calc_force_y + (node_mass * total_mass * 0.00000000006673 * (mass_of_y - node_coordinates_y) / pow(length_between_mass, 3));
        }
        else
        {

            if (tree->c0 != 0)
            { x_force_array[0], y_force_array[0] = getForce(node_n, tree->c0);  }
            if (tree->c1 != 0)
            { x_force_array[1], y_force_array[1] = getForce(node_n, tree->c1);  }
            if (tree->c2 != 0)
            { x_force_array[2], y_force_array[2] = getForce(node_n, tree->c2);  }
            if (tree->c3 != 0)
            { x_force_array[3], y_force_array[3] = getForce(node_n, tree->c3); }
           
	 	calc_force_x = x_force_array[0] + x_force_array[1] + x_force_array[2] + x_force_array[3];
		calc_force_y = y_force_array[0] + y_force_array[1] + y_force_array[2] + y_force_array[3];
        }
    }
    
	ret->x = calc_force_x;
	ret->y = calc_force_y;

return ret;

}

/*void simulat(struct quad_tree* q){

	struct node_list *hold_list = (struct node_list*) malloc(sizeof(struct node_list));
	struct quad_tree *hold_tree = init_quadTree(-536870911, 536870912, -536870911, 536870912);

	for (int i  = 0; i < nodesInTree->num_nodes; ++i){
		struct node *hold = nodesInTree->node_info[i]; 
		struct node *hold1 = (struct node*) malloc(sizeof(struct node));
		hold1->mass = hold->mass;
		double hold_xf;
		double hold_yf;
		struct pair* hold_f = get_force(q, hold);

		// compute force on node
		//hold_f = get_force(q, hold);
		hold_xf = hold_f->x;
		hold_yf = hold_f->y;
		free(hold_f);
		

		//get changes in velocity		
		hold1->x_velocity = hold_xf/hold1->mass;
		hold1->y_velocity = hold_yf/hold1->mass;

		//get changes in position
		hold1->x_coordinates = hold->x_coordinates+hold1->x_velocity;
		hold1->y_coordinates = hold->y_coordinates+hold1->y_velocity;
		hold1->x = hold1->x_coordinates*1000000;
		hold1->y = hold1->y_coordinates*1000000;

		//add to new list
		hold_list->node_info[i] = hold1;
		hold_list->num_nodes = i;
		// add to new tree
		insert(hold_tree, hold1);

		//if (i == 20 )
			//exit(0);
		
	}

	free(nodesInTree);
	free(q);
	nodesInTree = hold_list;
	q = hold_tree;
	set_up(q);
}

struct pair* get_force(struct quad_tree *q, struct node *n){
	double force_x = 0;
	double force_y = 0;
	struct pair *ret = (struct pair*) malloc(sizeof(struct pair));

	if(q->size == 1 || n->mass/sqrt(pow(q->cent_M_x - n->x_coordinates, 2) + pow(q->cent_M_y - n->y_coordinates, 2)) <= 1){
		force_x = 0.0000000000667*q->mass*(n->x_coordinates- q->cent_M_x)/pow(n->x_coordinates- q->cent_M_x, 3);
		force_y = 0.0000000000667*q->mass*(n->y_coordinates- q->cent_M_y)/pow(n->y_coordinates- q->cent_M_y, 3);
		ret->x = force_x;
		ret->y = force_y;
		return ret;
	}
	struct pair* hold;

	if(q->c0 != NULL){
		hold = get_force(q->c0,n);
		force_x += hold->x;
		force_y += hold->y;
	}
	if(q->c1 != NULL){
		hold = get_force(q->c1,n);
		force_x += hold->x;
		force_y += hold->y;
	}
	if(q->c2 != NULL){
		hold = get_force(q->c2,n);
		force_x += hold->x;
		force_y += hold->y;
	}
	if(q->c3 != NULL){
		hold = get_force(q->c3,n);
		force_x += hold->x;
		force_y += hold->y;
	}
	free(hold);

	ret->x = force_x;
	ret->y = force_y;
	return ret;
}*/
