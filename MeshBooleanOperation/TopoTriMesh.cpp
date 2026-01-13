#include "TopoTriMesh.h"
#include <set>

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

void TopoTriMesh::AddFace2TopoTriMesh(std::vector<Vec3d> const& pnts)
{
    // 1. 创建或获取顶点 v1, v2, v3
    Vertex* vertices[3] = { nullptr, nullptr, nullptr };
    for (size_t i = 0; i < pnts.size(); ++i) {
        auto it = p2V.find(pnts[i]);
        if (it != p2V.end()) {
            vertices[i] = it->second;
        }
        else {
            Vertex* newV = new Vertex{ pnts[i], nullptr, 0 };
            vs.push_back(newV);
            p2V[pnts[i]] = newV;
            vertices[i] = newV;
        }
    }

    // 如果无法获得3个顶点，返回
    if (!vertices[0] || !vertices[1] || !vertices[2]) return;

    // 2. 创建或获取边 e1, e2, e3
    auto makeEdgeKey = [](Vertex* a, Vertex* b) -> std::pair<Vertex*, Vertex*> {
        return std::make_pair(std::min(a, b), std::max(a, b));
        };

    Edge* edges[3] = { nullptr, nullptr, nullptr };
    std::map<std::pair<Vertex*, Vertex*>, Edge*> edgeMap;
    // 填充 edgeMap
    for (auto& e : es) {
        if (e && e->v1 && e->v2) {
            auto key = makeEdgeKey(e->v1, e->v2);
            edgeMap[key] = e;
        }
    }

    // 构造三条边：v1-v2, v2-v3, v3-v1
    for (int i = 0; i < 3; ++i) {
        auto key = makeEdgeKey(vertices[i], vertices[(i + 1) % 3]);
        if (edgeMap.count(key)) {
            edges[i] = edgeMap[key];
        }
        else {
            Edge* newE = new Edge{ vertices[i], vertices[(i + 1) % 3], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
            es.push_back(newE);
            edgeMap[key] = newE;
            edges[i] = newE;
        }
    }

    // 3. 设置边的拓扑关系（根据图片算法）
    Face* newF = new Face{ edges[0] };
    fs.push_back(newF);

    // 设置每条边的 lF, rF, lPE, rPE
    for (int i = 0; i < 3; ++i) {
        int next = (i + 1) % 3;
        int prev = (i - 1 + 3) % 3;

        if (!edges[i]->lF) {
            edges[i]->lF = newF;
            edges[i]->lPE = edges[prev];
            edges[i]->lSE = edges[next];
        }
        else {
            edges[i]->rF = newF;
            edges[i]->rPE = edges[prev];
            edges[i]->rSE = edges[next];
        }
    }

    // 确保每个顶点至少有一个关联的边
    for (int i = 0; i < 3; ++i) {
        if (!vertices[i]->e) {
            vertices[i]->e = edges[(i - 1 + 3) % 3];
        }
    }
}

void TopoTriMesh::RemoveFace(Face* f)
{
    auto RmEdgeFaceRel = [&](Edge* e) {
        if (!e) return;
        if (e->lF == f)
            e->lF = nullptr;
        else if (e->rF == f)
            e->rF = nullptr;
    };

    if (!f) return;
    Edge* fe = f->e, *pe = nullptr, *se = nullptr;
    if (!fe) return;
    if (fe->lF == f)
    {
        fe->lF = nullptr;
        pe = fe->lPE, se = fe->lSE;
    }
    else if (fe->rF == f)
    {
        fe->rF = nullptr;
        pe = fe->rPE, se = fe->rSE;
    }
    RmEdgeFaceRel(pe);
    RmEdgeFaceRel(se);

    auto newEnd = std::remove_if(fs.begin(), fs.end(), [&](Face* rmF) { return rmF == f; });
    fs.erase(newEnd, fs.end());
}

void TopoTriMesh::RemoveEdge(Edge* e)
{
    if (!e) return;

    auto cleanVertex = [&](Vertex* v) {
        if (!v) return;
        std::vector<Edge*> adjEs = v->GetAdjacentEdges();
        // 如果该顶点只连着这条边（即移除后孤立），则删除顶点
        if (adjEs.size() == 1 && adjEs[0] == e) {
            vs.erase(std::remove(vs.begin(), vs.end(), v), vs.end());
        }
        // 清理其他相邻边中指向 e 的拓扑指针
        for (Edge* adj : adjEs) {
            if (adj == e) continue; // 跳过自身
            if (adj->lPE == e) adj->lPE = nullptr;
            if (adj->lSE == e) adj->lSE = nullptr;
            if (adj->rPE == e) adj->rPE = nullptr;
            if (adj->rSE == e) adj->rSE = nullptr;
        }
        };

    cleanVertex(e->v1);
    cleanVertex(e->v2);

    // 从边列表中移除 e
    es.erase(std::remove(es.begin(), es.end(), e), es.end());
}

void TopoTriMesh::ReleaseMem()
{
    for (auto& f : fs)
    {
        delete f; 
        f = nullptr;
    }

    for (auto& e : es)
    {
        delete e;
        e = nullptr;
    }

    for (auto& v : vs)
    {
        delete v;
        v = nullptr;
    }
    p2V.clear();
}

std::vector<Edge*> Vertex::GetAdjacentEdges()
{
    if (!e) return {};

    std::vector<Edge*> cEdges;
    Edge* curEdge = e;

    // 保存初始边，防止无限循环
    Edge* startEdge = curEdge;

    do {
        cEdges.push_back(curEdge);

        // 根据 curEdge 是否以 v 为起点或终点，决定使用 lSE 或 rSE
        if (curEdge->v2 == this) {  // curEdge 的终点是 v → 使用 lSE
            curEdge = curEdge->lSE;
        }
        else if (curEdge->v1 == this) {  // curEdge 的起点是 v → 使用 rSE
            curEdge = curEdge->rSE;
        }
        else {
            // 理论上不会发生，但为了安全
            break;
        }

    } while (curEdge && curEdge != startEdge);

    return cEdges;
}

std::vector<Face*> Vertex::GetAdjacentFaces()
{
    std::set<Face*> faces;
    for (Edge* e : GetAdjacentEdges()) {
        if (e->lF) faces.insert(e->lF);
        if (e->rF) faces.insert(e->rF);
    }
    return std::vector<Face*>(faces.begin(), faces.end());
}

std::vector<Edge*> Face::getEdges()
{
    std::vector<Edge*> es(3);
    Edge* fe = e;
    es[1] = fe;
    if (!fe) return {};
    if (fe->lF == this)
    {
        if (fe->lPE)
            es[0] = fe->lPE;
        if (fe->lSE)
            es[2] = fe->lSE;
    }
    else
    {
        if (fe->rPE)
            es[0] = fe->rPE;
        if (fe->rSE)
            es[2] = fe->rSE;
    }
    return es;
}

std::vector<Vec3d> Face::getPnts()
{
    std::vector<Edge*> es = getEdges();
    std::vector<Vec3d> pnts;
    for (auto& e : es)
    {
        if (!e) continue;
        if (e->lF == this)
            pnts.push_back(e->v1->pnt);
        else
            pnts.push_back(e->v2->pnt);
    }
    return pnts;
}

std::vector<Vec3d> Edge::getPnts(bool left)
{
    if (!v1 || !v2) return {};
    if (left)
        return std::vector<Vec3d>{ v1->pnt, v2->pnt };
    return std::vector<Vec3d>{ v2->pnt, v1->pnt };
}