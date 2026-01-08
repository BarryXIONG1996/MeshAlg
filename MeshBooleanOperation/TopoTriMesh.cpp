#include "TopoTriMesh.h"

BndBox3d TopoTriMesh::GetBndBox() {
    BndBox3d bbox;
    if (vs.empty()) return bbox;

    bbox.lowerBnd = vs[0]->pnt;
    bbox.upperBnd = vs[0]->pnt;

    for (auto v : vs) {
        bbox.Add(v->pnt);
    }
    return bbox;
}