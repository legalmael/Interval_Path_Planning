#include "node.h"

Node::Node()
{
    distance=Infinity;
}

Node::~Node(){

}

void BisectNode(Node*& pX, Node*& pX1, Node*& pX2){

    pX1->space = box(2);  pX2->space = box(2);
    Bisect (pX->space,pX1->space,pX2->space);
    Node* pN(0);
    for(int it = 0; it != pX->Neighbors.size(); it++){
        pN = pX->Neighbors.at(it);
        for(int j = 0; j != pN->Neighbors.size(); j++){
            if (pN->Neighbors.at(j) == pX){
                pN->Neighbors.erase(pN->Neighbors.begin()+j);
                break;
            }
        }
        if (!Disjoint(pN->space, pX1->space)){ pX1->Neighbors.push_back(pN); pN->Neighbors.push_back(pX1);}
        if (!Disjoint(pN->space, pX2->space)){ pX2->Neighbors.push_back(pN); pN->Neighbors.push_back(pX2);}
    }
    pX1->Neighbors.push_back(pX2); pX2->Neighbors.push_back(pX1);
}
