#ifndef NODE_H
#define NODE_H

#include <iostream>
#include "vibes.h"
#include "interval.h"
#include "box.h"
#include <list>

#define Infinity 10000

enum iboolean{itrue, ifalse, iperhaps};

class Node
{
public:
    Node();
    ~Node();

    box space;
    iboolean category;      // Classification of the configuration space
    vector<Node*> Neighbors;   // Graph representation
    int distance;           // Shortest path algorithm

};

    void BisectNode(Node*& pX, Node*& pX1, Node*& pX2);

#endif // NODE_H
