#include "BooleanOpHelper.h"
#include <queue>
#include <set>

// 构造函数
BooleanOpHelper::BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, TopoTriMesh& coPlanes, int opType)
    : m_objMesh(objMesh), m_subMesh(subMesh), m_coPlanes(coPlanes), m_opType(opType) {}

// 主执行函数
bool BooleanOpHelper::Execute(TopoTriMesh& res) {
    // Step 1: 提取区域
    TopoTriMesh Mooutt, Moint;
    BFSExtractRegion(m_objMesh, Mooutt, Moint);

    TopoTriMesh Mtouto, Mtino;
    BFSExtractRegion(m_subMesh, Mtouto, Mtino);

    // Step 2: 共面部分
    TopoTriMesh Mtono = m_coPlanes;

    // Step 3: 按操作类型组合
    std::vector<TopoTriMesh> ms;
    switch (m_opType) {
    case 0: // INTERSECTION
        res = Moint;
        ms = { Mtino, Mtono };
        ReleaseMeshExceptBoundary(Mooutt);
        ReleaseMeshExceptBoundary(Mtouto);
        CombineTopoTriMesh(res, ms);
        break;
    case 1: // UNION
        res = Mooutt;
        ms = { Mtouto, Mtono };
        ReleaseMeshExceptBoundary(Moint);
        ReleaseMeshExceptBoundary(Mtino);
        CombineTopoTriMesh(res, ms);
        break;
    case 2: // DIFFERENCE (obj - sub)
        res = Mooutt;
        ms = { Mtino };
        ReleaseMeshExceptBoundary(Moint);
        ReleaseMeshExceptBoundary(Mtouto);
        Mtono.ReleaseMem();
        CombineTopoTriMesh(res, ms);
        break;
    default:
        return false;
    }
    return true;
}

// BFS 区域提取
void BooleanOpHelper::BFSExtractRegion(const TopoTriMesh& mesh, TopoTriMesh& Mout, TopoTriMesh& Min) {
    std::queue<Edge*> es;
    std::set<Edge*> visitedE;
    std::set<Vertex*> vvOut;
    std::set<Edge*> veOut;
    std::set<Face*> vfOut;

    // 初始化：从 posTag == 2 (Out) 的顶点出发
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
    if (!e || veOut.find(e) != veOut.end()) return;

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

void BooleanOpHelper::CombineTopoTriMesh(TopoTriMesh& M, std::vector<TopoTriMesh>& Ms)
{
    for (auto& m : Ms)
    {
        for (auto& f : m.fs) // 遍历M.fs(f) :
        {
            if (!f) continue;
            Edge* es[3];
            Edge* fe = f->e;
            es[1] = fe;
            if (!fe) continue;
            if (fe->lF == f)
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
            std::vector<Vec3d> pnts;
            for (auto& e : es)
            {
                if (!e) continue;
                if (e->lF == f)
                    pnts.push_back(e->v1->pnt);
                else
                    pnts.push_back(e->v2->pnt);
            }
            if (pnts.size() < 3) continue;
            // 将pnts添加到M
            M.AddFace2TopoTriMesh(pnts);
        }
        // 内存释放
        std::for_each(m.fs.begin(), m.fs.end(), [](Face* f) { delete f; f = nullptr; });
        std::for_each(m.es.begin(), m.es.end(), [](Edge* e) { delete e; e = nullptr; });
        std::for_each(m.vs.begin(), m.vs.end(), [](Vertex* v) { delete v; v = nullptr; });
    }
}

void BooleanOpHelper::ReleaseMeshExceptBoundary(TopoTriMesh& M)
{
    // 遍历M
    std::set<Face*> fs;
    std::set<Edge*> es, esOn;
    for (auto& f : M.fs)
    {
        if (!f) continue;
        fs.insert(f);
    }

    for (auto& e : M.es)
    {
        if (!e) continue;
        if (e->v1 && e->v2 && e->v1->posTag == 3 && e->v2->posTag == 3)
            esOn.insert(e);
        else
            es.insert(e);
    }

    // 解除顶点联系
    for (auto& v : M.vs)
    {
        if (!v || v->posTag != 3 || !es.count(v->e)) continue;
        // 获取与v相连的所有边
        std::vector<Edge*> cEdges = v->GetAdjacentEdges();
        for (auto& ce : cEdges)
        {
            if (ce->v2 == v && !es.count(ce))
            {
                v->e = ce;
                break;
            }
        }
    }

    // 解除交线联系
    for (auto& e : esOn)
    {
        if (fs.count(e->lF))
        {
            e->lF = nullptr;
            e->lSE = e->lPE = nullptr;
        }
        else
        {
            e->rF = nullptr;
            e->rSE = e->rPE = nullptr;
        }
    }

    // 内存释放
    std::for_each(M.fs.begin(), M.fs.end(), [](Face* f) { delete f; f = nullptr; });
    std::for_each(M.es.begin(), M.es.end(), [](Edge* e) {
        if (e->v1 && e->v2 && e->v1->posTag == 3 && e->v2->posTag == 3)
            return;
        delete e;
        e = nullptr;
        });
    std::for_each(M.vs.begin(), M.vs.end(), [](Vertex* v) {
        if (v->posTag == 3)
            return;
        delete v;
        v = nullptr;
        });
}
