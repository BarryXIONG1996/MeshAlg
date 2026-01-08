#include "BooleanOpHelper.h"
#include <queue>
#include <set>

// 构造函数
BooleanOpHelper::BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, TopoTriMesh& coPlanes, int opType)
    : m_objMesh(objMesh), m_subMesh(subMesh), m_coPlanes(coPlanes), m_opType(opType) {}

// 主执行函数
bool BooleanOpHelper::Execute(TopoTriMesh& res) {
    // Step 1: 提取区域（变量名按意见(4)修正）
    TopoTriMesh Mooutt, Moint;
    BFSExtractRegion(m_objMesh, Mooutt, Moint);

    TopoTriMesh Mtouto, Mtino;
    BFSExtractRegion(m_subMesh, Mtouto, Mtino);

    // Step 2: 共面部分
    TopoTriMesh Mtono = m_coPlanes;

    // Step 3: 按操作类型组合
    switch (m_opType) {
    case 0: // INTERSECTION
        res = CombineMeshes(Moint, Mtino, Mtono);
        break;
    case 1: // UNION
        res = CombineMeshes(Mooutt, Mtouto, Mtono);
        break;
    case 2: // DIFFERENCE (obj - sub)
        res = CombineMeshes(Mooutt, Mtino, TopoTriMesh{});
        break;
    default:
        return false;
    }
    return true;
}

// BFS 区域提取（核心修正点集中于此）
void BooleanOpHelper::BFSExtractRegion(const TopoTriMesh& mesh, TopoTriMesh& Mout, TopoTriMesh& Min) {
    std::queue<Edge*> es;
    std::set<Edge*> visitedE; // ← 修改(1): 改为 set
    std::set<Vertex*> vvOut;
    std::set<Edge*> veOut;
    std::set<Face*> vfOut;

    // 初始化：从 posTag == 2 (Out) 的顶点出发 ← 修改(3)
    for (auto* v : mesh.vs) {
        if (v && v->posTag == 2) {
            Edge* e = v->e;
            if (e && visitedE.find(e) == visitedE.end()) {
                es.push(e);
                visitedE.insert(e);
            }
        }
    }

    // BFS 主循环
    while (!es.empty()) {
        Edge* e = es.front();
        es.pop();

        AddFace(e->lF, Mout, vvOut, veOut, vfOut);
        AddFace(e->rF, Mout, vvOut, veOut, vfOut);

        AddAdjacentEdges(e->lF, e, es, visitedE);
        AddAdjacentEdges(e->rF, e, es, visitedE);
    }

    // 构建 Min：未被 Mout 收录的部分
    for (auto* f : mesh.fs) {
        if (vfOut.find(f) == vfOut.end()) {
            Min.fs.push_back(f);
        }
    }

    // ← 修改(2): 移除 Min.es 的 std::find 重复检查，直接 push
    for (auto* e : mesh.es) {
        if (veOut.find(e) == veOut.end() ||
            (e->v1 && e->v2 && e->v1->posTag == 3 && e->v2->posTag == 3)) {
            Min.es.push_back(e); // 允许少量重复，或由上层处理
        }
    }

    for (auto* v : mesh.vs) {
        // ON 顶点（posTag==3）强制加入 Min ← 修改(3)
        if (vvOut.find(v) == vvOut.end() || v->posTag == 3) {
            Min.vs.push_back(v);
        }
    }
}

// AddFace：保持原始三边数组逻辑，无修改
void BooleanOpHelper::AddFace(Face* f, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve, std::set<Face*>& vf) {
    if (!f || vf.find(f) != vf.end()) return;

    M.fs.push_back(f);
    vf.insert(f);

    Edge* fe = f->e;
    AddEdge(fe, M, vv, ve);

    if (fe->lF == f) {
        AddEdge(fe->lPE, M, vv, ve);
        AddEdge(fe->lSE, M, vv, ve);
    }
    else {
        AddEdge(fe->rPE, M, vv, ve);
        AddEdge(fe->rSE, M, vv, ve);
    }
}

// AddEdge：保持原样
void BooleanOpHelper::AddEdge(Edge* e, TopoTriMesh& Mout, std::set<Vertex*>& vvOut, std::set<Edge*>& veOut) {
    if (!e) return;
    if (veOut.find(e) != veOut.end()) return;

    Mout.es.push_back(e);
    veOut.insert(e);

    if (e->v1 && vvOut.find(e->v1) == vvOut.end()) {
        Mout.vs.push_back(e->v1);
        vvOut.insert(e->v1);
    }
    if (e->v2 && vvOut.find(e->v2) == vvOut.end()) {
        Mout.vs.push_back(e->v2);
        vvOut.insert(e->v2);
    }
}

// AddAdjacentEdges：仅增加 ON 边跳过条件 ← 修改(5)
void BooleanOpHelper::AddAdjacentEdges(Face* f, Edge* skipEdge, std::queue<Edge*>& es, std::set<Edge*>& visitedE) {
    if (!f) return;

    Edge* fe = f->e;
    Edge* edgesToCheck[3];
    edgesToCheck[0] = fe;
    if (fe->lF == f) {
        edgesToCheck[1] = fe->lPE;
        edgesToCheck[2] = fe->lSE;
    }
    else {
        edgesToCheck[1] = fe->rPE;
        edgesToCheck[2] = fe->rSE;
    }

    for (Edge* adjE : edgesToCheck) {
        // ← 修改(5): 跳过 ON 边（两个端点 posTag == 3）
        if (adjE && adjE->v1 && adjE->v2 &&
            adjE->v1->posTag == 3 && adjE->v2->posTag == 3) {
            continue;
        }

        if (adjE && adjE != skipEdge && !visitedE.count(adjE)) {
            es.push(adjE);
            visitedE.insert(adjE);
        }
    }
}

// CombineMeshes：实现实际的合并逻辑
TopoTriMesh BooleanOpHelper::CombineMeshes(const TopoTriMesh& mesh1, const TopoTriMesh& mesh2, const TopoTriMesh& mesh3) {
    TopoTriMesh result;

    // ToDo: 实现实际的合并逻辑

    return result;
}