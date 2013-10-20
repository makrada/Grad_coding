/*
 * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
 *
 ******         This file contains the basic structs of the program  **********
 **/

struct node{
  struct node *nxt;
  struct elem *element;
};

struct dc_nodes{
  char  *type;
  double start;
  double end;
  double step;
  int pos;
  int pos2;
  char  *tplot;
  int value_plot;
  double *trans_spec;
  int flag;

};

struct node1{
  struct node1 *nxt;
  struct dc_nodes *element;
};

struct elem{
  char *type;
  double *term;
  double *value;
};
