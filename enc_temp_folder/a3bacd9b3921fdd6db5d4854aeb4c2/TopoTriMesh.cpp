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

Face* TopoTriMesh::AddFace2TopoTriMesh(std::vector<Vec3d> const& pnts)
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
    if (!vertices[0] || !vertices[1] || !vertices[2]) return nullptr;

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
    newF->topo = this;
    fs.push_back(newF);
    for (auto const& pnt : pnts) newF->bbox.Add(pnt);

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

    return newF;
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
        if (v->e == e) {
            for (Edge* adj : adjEs) {
                if (adj == e || adj->v2 != v)
                    continue;
                v->e = adj;
                break;
            }
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
    fs.clear();

    for (auto& e : es)
    {
        delete e;
        e = nullptr;
    }
    es.clear();

    for (auto& v : vs)
    {
        delete v;
        v = nullptr;
    }
    es.clear();

    p2V.clear();
}

void TopoTriMesh::Build(const TriMesh& triMesh) {
    ReleaseMem();

    if (triMesh.indices.empty() || triMesh.points.empty()) {
        return;
    }

    // --- Step 1: 构建去重顶点 ---
    std::map<Vertex*, int> v2Idx;           // Vertex* -> index in vs
    std::map<int, int> pIdx2vIdx;           // original point index -> vs index

    for (int idx = 0; idx < static_cast<int>(triMesh.points.size()); ++idx) {
        const Vec3d& pnt = triMesh.points[idx];
        auto it = p2V.find(pnt);
        if (it == p2V.end()) { // 新顶点
            Vertex* v = new Vertex{ pnt };
            int vidx = static_cast<int>(vs.size());
            vs.push_back(v);
            p2V[pnt] = v;
            v2Idx[v] = vidx;
            pIdx2vIdx[idx] = vidx;
        }
        else { // 已存在
            Vertex* v = it->second;
            int vidx = v2Idx[v];
            pIdx2vIdx[idx] = vidx;
        }
    }

    // --- Step 2: 构建面和边 ---
    std::map<std::pair<int, int>, Edge*> cEs; // (min, max) -> Edge*

    for (int idx = 0; idx + 3 <= static_cast<int>(triMesh.indices.size()); idx += 4) {
        // 跳过无效三角形（如索引为0）
        if (triMesh.indices[idx] == 0 ||
            triMesh.indices[idx + 1] == 0 ||
            triMesh.indices[idx + 2] == 0) {
            continue;
        }

        Face* f = new Face;
        fs.push_back(f);
        f->topo = this;

        // 1-based → 0-based point index → vs index
        int i0 = pIdx2vIdx.at(triMesh.indices[idx] - 1);
        int i1 = pIdx2vIdx.at(triMesh.indices[idx + 1] - 1);
        int i2 = pIdx2vIdx.at(triMesh.indices[idx + 2] - 1);

        f->bbox.Add(vs.at(i0)->pnt);
        f->bbox.Add(vs.at(i1)->pnt);
        f->bbox.Add(vs.at(i2)->pnt);

        // 构造边 key
        std::pair<int, int> e1_key = { std::min(i0, i1), std::max(i0, i1) };
        std::pair<int, int> e2_key = { std::min(i1, i2), std::max(i1, i2) };
        std::pair<int, int> e3_key = { std::min(i2, i0), std::max(i2, i0) };

        Edge* tes[3];
        bool isENew[3] = { false, false, false };

        // --- 处理 e0: (i0, i1) ---
        auto it1 = cEs.find(e1_key);
        if (it1 == cEs.end()) {
            tes[0] = new Edge{ vs[i0], vs[i1], f };
            isENew[0] = true;
            es.push_back(tes[0]);
            cEs[e1_key] = tes[0];
            if (!vs[i1]->e) vs[i1]->e = tes[0];
        }
        else {
            tes[0] = it1->second;
            tes[0]->rF = f; // 第二个邻接面
        }
        f->e = tes[0];

        // --- 处理 e1: (i1, i2) ---
        auto it2 = cEs.find(e2_key);
        if (it2 == cEs.end()) {
            tes[1] = new Edge{ vs[i1], vs[i2], f };
            isENew[1] = true;
            es.push_back(tes[1]);
            cEs[e2_key] = tes[1];
            if (!vs[i2]->e) vs[i2]->e = tes[1];
        }
        else {
            tes[1] = it2->second;
            tes[1]->rF = f;
        }

        // --- 处理 e2: (i2, i0) ---
        auto it3 = cEs.find(e3_key);
        if (it3 == cEs.end()) {
            tes[2] = new Edge{ vs[i2], vs[i0], f };
            isENew[2] = true;
            es.push_back(tes[2]);
            cEs[e3_key] = tes[2];
            if (!vs[i0]->e) vs[i0]->e = tes[2];
        }
        else {
            tes[2] = it3->second;
            tes[2]->rF = f;
        }

        // --- 设置边的邻接边（绕面顺序）---
        // 对于新边：设置 lPE (left previous edge), lSE (left next edge)
        // 对于已有边：设置 rPE (right previous), rSE (right next)
        for (int i = 0; i < 3; ++i) {
            Edge* cur = tes[i];
            Edge* prev = tes[(i - 1 + 3) % 3]; // 前一条边（逆时针）
            Edge* next = tes[(i + 1) % 3];     // 后一条边

            if (isENew[i]) {
                // 当前面在边的 "左侧"
                cur->lPE = prev;
                cur->lSE = next;
            }
            else {
                // 当前面在边的 "右侧"
                cur->rPE = prev;
                cur->rSE = next;
            }
        }
    }
}

void TopoTriMesh::ToMesh(TriMesh& mesh)
{
    mesh.points.clear();
    mesh.indices.clear();

    if (vs.empty()) {
        return;
    }

    // Step 1: 导出所有顶点，并建立 Vertex* → 1-based index 映射
    std::map<Vertex*, int> vertexToIndex;
    mesh.points.reserve(vs.size());
    for (size_t i = 0; i < vs.size(); ++i) {
        mesh.points.push_back(vs[i]->pnt);
        vertexToIndex[vs[i]] = static_cast<int>(i + 1); // 1-based
    }

    // Step 2: 遍历每个面，导出三角形（带 0 分隔符）
    for (Face* f : fs) {
        if (!f) continue;

        std::vector<Vertex*> verts = f->getVertices();
        if (verts.size() != 3) continue; // 只处理三角形

        // 获取 1-based 索引
        auto it0 = vertexToIndex.find(verts[0]);
        auto it1 = vertexToIndex.find(verts[1]);
        auto it2 = vertexToIndex.find(verts[2]);
        if (it0 == vertexToIndex.end() ||
            it1 == vertexToIndex.end() ||
            it2 == vertexToIndex.end()) {
            continue; // 安全检查
        }

        // 写入：i0, i1, i2, 0
        mesh.indices.push_back(it0->second);
        mesh.indices.push_back(it1->second);
        mesh.indices.push_back(it2->second);
        mesh.indices.push_back(0);
    }
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

std::vector<Vertex*> Face::getVertices()
{
    std::vector<Edge*> es = getEdges();
    std::vector<Vertex*> pnts;
    for (auto& e : es)
    {
        if (!e) continue;
        if (e->lF == this)
            pnts.push_back(e->v1);
        else
            pnts.push_back(e->v2);
    }
    return pnts;
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