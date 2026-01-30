#include "BooleanOpHelper.h"
#include <queue>
#include <set>

extern const double g_epsilon;

// 构造函数
BooleanOpHelper::BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, TopoTriMesh& coPlanes, int opType)
    : m_objMesh(objMesh), m_subMesh(subMesh), m_coPlanes(coPlanes), m_opType(opType) {}

void ShowTriMesh(const TriMesh& mesh);

// 主执行函数
bool BooleanOpHelper::Execute(TopoTriMesh& res) {
    // Step 1: 提取区域
    TopoTriMesh Mooutt, Moint;
    BFSExtractRegion(m_objMesh, Mooutt, Moint);

    TopoTriMesh Mtouto, Mtino;
    BFSExtractRegion(m_subMesh, Mtouto, Mtino);
    
#ifdef _DEBUG
    //TriMesh oout, oint, touto, tino;
    TriMesh oot, tio;
    //Mooutt.ToMesh(oout);
    //Moint.ToMesh(oint);
    //Mtouto.ToMesh(touto);
    //Mtino.ToMesh(tino);
    //ShowTriMesh(oout);
    //ShowTriMesh(oint);
    //ShowTriMesh(touto);
    //ShowTriMesh(tino);
#endif

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
        Mtino.Reverse();
        ms = { Mtino };
        ReleaseMeshExceptBoundary(Moint);
        ReleaseMeshExceptBoundary(Mtouto);
        Mtono.ReleaseMem();


#ifdef _DEBUG
        Mooutt.ToMesh(oot);
        Mtino.ToMesh(tio);
        ShowTriMesh(oot);
        ShowTriMesh(tio);
#endif

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
            for (Edge* e : v->es) {
                if (e && visitedE.find(e) == visitedE.end()) {
                    es.push(e);
                    visitedE.insert(e);
                }
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

#if 1 // 避免出现悬空边和顶点
    std::set<Vertex*> vvIn;
    std::set<Edge*> veIn;
    for (auto* f : Min.fs) {
        std::vector<Edge*> fes = f->getEdges();
        std::vector<Vertex*> fvs = f->getVertices();
        for (Edge* e : fes) {
            if (veIn.count(e)) continue;
            Min.es.push_back(e);
            veIn.insert(e);
        }
        for (Vertex* v : fvs) {
            if (vvIn.count(v)) continue;
            Min.vs.push_back(v);
            vvIn.insert(v);
        }
    }
#else
    /*目前的逻辑可能导致生成的拓扑实体中包含悬空边（即不与面关联的边）*/
    for (auto* e : mesh.es) {
        if (veOut.find(e) == veOut.end() ||
            (e->v1 && e->v2 && e->v1->posTag == 3 && e->v2->posTag == 3)) {
            Min.es.push_back(e); // 允许少量重复，或由上层处理
        }
    }

    /*目前的逻辑可能导致生成的拓扑实体中包含悬空顶点（即不与面关联的顶点）*/
    for (auto* v : mesh.vs) {
        // ON 顶点（posTag==3）强制加入 Min ← 修改(3)
        if (vvOut.find(v) == vvOut.end() || v->posTag == 3) {
            Min.vs.push_back(v);
            Min.p2V[v->pnt] = v;
        }
    }
#endif
}

// AddFace：保持原始三边数组逻辑，无修改
void BooleanOpHelper::AddFace(Face* f, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve, std::set<Face*>& vf) {
    if (!f || vf.find(f) != vf.end()) return;

    M.fs.push_back(f);
    vf.insert(f);

    Edge* fe = f->es[0];
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
        Mout.p2V[e->v1->pnt] = e->v1;
        vvOut.insert(e->v1);
    }
    if (e->v2 && vvOut.find(e->v2) == vvOut.end()) {
        Mout.vs.push_back(e->v2);
        Mout.p2V[e->v2->pnt] = e->v2;
        vvOut.insert(e->v2);
    }
}

// AddAdjacentEdges：仅增加 ON 边跳过条件 ← 修改(5)
void BooleanOpHelper::AddAdjacentEdges(Face* f, Edge* skipEdge, std::queue<Edge*>& es, std::set<Edge*>& visitedE) {
    if (!f) return;

    Edge* fe = f->es[0];
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
            std::vector<Vec3d> pnts = f->getPnts();
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

void BooleanOpHelper::ReleaseMeshExceptBoundary(TopoTriMesh& M) {
    // 遍历M
    std::set<Face*> fs;
    for (auto& f : M.fs) {
        if (!f) continue;
        fs.insert(f);
    }

    // 是否是内部面
    auto IsInnerEdges = [&](Edge* e) {
        if (e->lF && !fs.count(e->lF)) return false;
        if (e->rF && !fs.count(e->rF)) return false;
        return true;
    };

    std::set<Edge*> es/*内部边*/, esOn/*边界边*/;
    for (auto& e : M.es) {
        if (!e) continue;
        if ((e->v1 && e->v2) 
            && (e->v1->posTag == 3 && e->v2->posTag == 3)
            && !IsInnerEdges(e)/*非内部边*/)
            esOn.insert(e);
        else
            es.insert(e);
    }

    std::set<Vertex*> vsOn/*边界顶点*/;
    for (auto& e : esOn) {
        // 解除边界联系
        if (fs.count(e->lF)) {
            e->lF = nullptr;
            e->lSE = e->lPE = nullptr;
        }
        else {
            e->rF = nullptr;
            e->rSE = e->rPE = nullptr;
        }
        // 解除边界顶点联系
        if (e->v1) {
            vsOn.insert(e->v1);
            std::vector<Edge*> cEdges = e->v1->GetAdjacentEdges();
            for (auto& ce : cEdges) {
                if (es.count(ce))
                    e->v1->es.erase(ce);
            }
        }
        if (e->v2) {
            vsOn.insert(e->v2);
            std::vector<Edge*> cEdges = e->v2->GetAdjacentEdges();
            for (auto& ce : cEdges) {
                if (es.count(ce))
                    e->v2->es.erase(ce);
            }
        }
    }

    // 内存释放
    std::for_each(M.fs.begin(), M.fs.end(), [](Face* f) { delete f; f = nullptr; });
    std::for_each(M.es.begin(), M.es.end(), [&](Edge* e) {
        if (esOn.count(e)) return; // 跳过边界
        delete e;
        e = nullptr;
    });
    std::for_each(M.vs.begin(), M.vs.end(), [&](Vertex* v) {
        if (vsOn.count(v)) return; // 跳过边界顶点
        delete v;
        v = nullptr;
    });
}
