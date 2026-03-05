#include <string>
struct Node {
  struct Node *left;
  struct Node *right;
  struct Node *parent;
  double distance;
  unsigned int nodenumber;
  unsigned int *child_ids;

  Node(Node *left, Node *right, Node *parent): left(left), right(right), parent(parent) {}
};

struct Tree {
  struct Tnode *child[3];
  unsigned int numnodes;
};

struct Cluster {
  unsigned int numclusters;
  struct Tnode *node;
  struct DistanceMatrix *matrix;
};
