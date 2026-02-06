#include <fstream>
#include <sstream>
#include <cctype>
#include <array>
#include <queue>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>
#include <map>

#include "Geometry.h"

double g_epsilon = 1e-12;

// ========== 构造函数定义 ==========
Mat4d::Mat4d() {
    SetIdentity();
}

Mat4d::Mat4d(const double* data) {
    std::copy(data, data + 16, m);
}

Mat4d::Mat4d(double m00, double m01, double m02, double m03,
    double m10, double m11, double m12, double m13,
    double m20, double m21, double m22, double m23,
    double m30, double m31, double m32, double m33)
{
    m[0] = m00; m[4] = m01; m[8] = m02; m[12] = m03;
    m[1] = m10; m[5] = m11; m[9] = m12; m[13] = m13;
    m[2] = m20; m[6] = m21; m[10] = m22; m[14] = m23;
    m[3] = m30; m[7] = m31; m[11] = m32; m[15] = m33;
}

// ========== 静态工厂函数定义 ==========
Mat4d Mat4d::Identity() {
    return Mat4d(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
}

Mat4d Mat4d::Translation(const Vec3d& t) {
    return Mat4d(1, 0, 0, t.x, 0, 1, 0, t.y, 0, 0, 1, t.z, 0, 0, 0, 1);
}

Mat4d Mat4d::Scaling(const Vec3d& s) {
    return Mat4d(s.x, 0, 0, 0, 0, s.y, 0, 0, 0, 0, s.z, 0, 0, 0, 0, 1);
}

Mat4d Mat4d::Rotation(const Vec3d& axis, double angle) {
    // 罗德里格斯旋转公式
    // 将向量分解为平行 / 垂直于旋转轴的两个分量
    // 旋转垂直分量（平面内旋转）
    Vec3d k = axis.Normalization();
    double c = std::cos(angle);
    double s = std::sin(angle);
    double omc = 1.0 - c;

    double xx = k.x * k.x, yy = k.y * k.y, zz = k.z * k.z;
    double xy = k.x * k.y, xz = k.x * k.z, yz = k.y * k.z;

    return Mat4d(
        xx * omc + c, xy * omc - k.z * s, xz * omc + k.y * s, 0,
        xy * omc + k.z * s, yy * omc + c, yz * omc - k.x * s, 0,
        xz * omc - k.y * s, yz * omc + k.x * s, zz * omc + c, 0,
        0, 0, 0, 1
    );
}

Mat4d Mat4d::RotationX(double angle) {
    /*y' = y·cosθ - z·sinθ
      z' = y·sinθ + z·cosθ*/
    double c = std::cos(angle), s = std::sin(angle);
    return Mat4d(1, 0, 0, 0, 0, c, -s, 0, 0, s, c, 0, 0, 0, 0, 1);
}

Mat4d Mat4d::RotationY(double angle) {
    double c = std::cos(angle), s = std::sin(angle);
    return Mat4d(c, 0, s, 0, 0, 1, 0, 0, -s, 0, c, 0, 0, 0, 0, 1);
}

Mat4d Mat4d::RotationZ(double angle) {
    double c = std::cos(angle), s = std::sin(angle);
    return Mat4d(c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
}

// ========== 成员函数定义 ==========
void Mat4d::SetIdentity() {
    std::fill(std::begin(m), std::end(m), 0.0);
    m[0] = m[5] = m[10] = m[15] = 1.0;
}

void Mat4d::Transpose() {
    std::swap(m[1], m[4]); std::swap(m[2], m[8]); std::swap(m[3], m[12]);
    std::swap(m[6], m[9]); std::swap(m[7], m[13]); std::swap(m[11], m[14]);
}

Mat4d Mat4d::operator*(const Mat4d& other) const {
    Mat4d res;
    for (int col = 0; col < 4; ++col) {
        double c0 = m[col * 4 + 0], c1 = m[col * 4 + 1], c2 = m[col * 4 + 2], c3 = m[col * 4 + 3];
        res.m[col * 4 + 0] = c0 * other.m[0] + c1 * other.m[1] + c2 * other.m[2] + c3 * other.m[3];
        res.m[col * 4 + 1] = c0 * other.m[4] + c1 * other.m[5] + c2 * other.m[6] + c3 * other.m[7];
        res.m[col * 4 + 2] = c0 * other.m[8] + c1 * other.m[9] + c2 * other.m[10] + c3 * other.m[11];
        res.m[col * 4 + 3] = c0 * other.m[12] + c1 * other.m[13] + c2 * other.m[14] + c3 * other.m[15];
    }
    return res;
}

Vec3d Mat4d::TransformPoint(const Vec3d& v) const {
    double x = m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12];
    double y = m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13];
    double z = m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14];
    double w = m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15];
    if (std::abs(w) > 1e-12) {
        double invW = 1.0 / w;
        x *= invW; y *= invW; z *= invW;
    }
    return Vec3d{ x, y, z };
}

Vec3d Mat4d::TransformVector(const Vec3d& v) const {
    return Vec3d{
        m[0] * v.x + m[4] * v.y + m[8] * v.z,
        m[1] * v.x + m[5] * v.y + m[9] * v.z,
        m[2] * v.x + m[6] * v.y + m[10] * v.z
    };
}

double& Mat4d::operator()(int row, int col) {
    assert(row >= 0 && row < 4 && col >= 0 && col < 4);
    return m[col * 4 + row];
}

double Mat4d::operator()(int row, int col) const {
    assert(row >= 0 && row < 4 && col >= 0 && col < 4);
    return m[col * 4 + row];
}

Vec2d Vec2d::operator+(const Vec2d& v) const {
    return { x + v.x, y + v.y };
}

Vec2d Vec2d::operator-(const Vec2d& v) const {
    return { x - v.x, y - v.y };
}

Vec2d Vec2d::operator*(double s) const {
    return { x * s, y * s };
}

Vec2d Vec2d::operator/(double s) const {
    return { x / s, y / s };
}

// 点积
double Vec2d::Dot(const Vec2d& v) const {
    return x * v.x + y * v.y;
}

// 2D 叉积（标量）
double Vec2d::Cross(const Vec2d& v) const {
    return x * v.y - y * v.x;
}

// 长度
double Vec2d::Length() const {
    return sqrt(x * x + y * y);
}

// 长度平方
double Vec2d::LengthSq() const {
    return x * x + y * y;
}

// 归一化
Vec2d Vec2d::Normalization() const {
    double len = Length();
    if (len  < g_epsilon) {
        return { 0.0, 0.0 };
    }
    return { x / len, y / len };
}

// 平行判断
bool Vec2d::Parallel(const Vec2d& v) const {
    return std::abs(Cross(v)) < g_epsilon;
}

// 相等判断
bool Vec2d::Equal(const Vec2d& v) const {
    return operator-(v).Length() < g_epsilon;
}

Vec3d Vec3d::operator-(const Vec3d& v) const {
    return Vec3d{ x - v.x, y - v.y, z - v.z };
}

Vec3d Vec3d::operator+(const Vec3d& v) const {
    return Vec3d{ x + v.x, y + v.y, z + v.z };
}

Vec3d Vec3d::operator*(double s) const {
    return Vec3d{ x * s, y * s, z * s };
}

Vec3d Vec3d::operator/(double s) const {
    return { x / s, y / s, z / s };
}

double Vec3d::Dot(const Vec3d& v) const {
    return x * v.x + y * v.y + z * v.z;
}

Vec3d Vec3d::Cross(const Vec3d& v) const {
    return Vec3d{
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x
    };
}

double Vec3d::Length() const {
    return sqrt(x * x + y * y + z * z);
}

double Vec3d::LengthSq() const
{
    return x * x + y * y + z * z;
}

Vec3d Vec3d::Normalization() const
{
    double len = Length();

    if (len < g_epsilon)
        return { 0,0,0 };

    return { x/len, y/len, z/len };
}

bool Vec3d::Parallel(const Vec3d& v) const
{
    Vec3d crossVec = Cross(v);
    return crossVec.Length() < g_epsilon;
}

bool Vec3d::Equal(const Vec3d& v) const
{
    return operator-(v).Length() < g_epsilon;
}

BndBox3d::BndBox3d()
{
    lowerBnd = { 
        std::numeric_limits<double>::max(), 
        std::numeric_limits<double>::max(), 
        std::numeric_limits<double>::max() 
    };

    upperBnd = {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };
}

BndBox3d BndBox3d::Intersect(BndBox3d const& bnd) const {
    BndBox3d res;
    res.lowerBnd.x = std::max(lowerBnd.x, bnd.lowerBnd.x);
    res.lowerBnd.y = std::max(lowerBnd.y, bnd.lowerBnd.y);
    res.lowerBnd.z = std::max(lowerBnd.z, bnd.lowerBnd.z);
    res.upperBnd.x = std::min(upperBnd.x, bnd.upperBnd.x);
    res.upperBnd.y = std::min(upperBnd.y, bnd.upperBnd.y);
    res.upperBnd.z = std::min(upperBnd.z, bnd.upperBnd.z);
    return res;
}

bool BndBox3d::IsOut(BndBox3d bnd) const {
    // 判断两个包围盒是否不相交(精度g_epsilon)
    return (upperBnd.x < bnd.lowerBnd.x - g_epsilon || lowerBnd.x > bnd.upperBnd.x + g_epsilon ||
        upperBnd.y < bnd.lowerBnd.y - g_epsilon || lowerBnd.y > bnd.upperBnd.y + g_epsilon ||
        upperBnd.z < bnd.lowerBnd.z - g_epsilon || lowerBnd.z > bnd.upperBnd.z + g_epsilon);
}

void BndBox3d::Add(const Vec3d& p) {
    if (p.x < lowerBnd.x) lowerBnd.x = p.x;
    if (p.y < lowerBnd.y) lowerBnd.y = p.y;
    if (p.z < lowerBnd.z) lowerBnd.z = p.z;
    if (p.x > upperBnd.x) upperBnd.x = p.x;
    if (p.y > upperBnd.y) upperBnd.y = p.y;
    if (p.z > upperBnd.z) upperBnd.z = p.z;
}

Vec3d BndBox3d::Center() const {
    return {
        (lowerBnd.x + upperBnd.x) / 2,
        (lowerBnd.y + upperBnd.y) / 2,
        (lowerBnd.z + upperBnd.z) / 2,
    };
}

void TriMesh::Transform(const Mat4d& mat)
{
    for (Vec3d& p : points) {
        p = mat.TransformPoint(p);
    }
}

void TriMesh::BuildFromOBJ(const std::string& fileName)
{
    // 清空现有数据
    points.clear();
    indices.clear();

    std::ifstream file(fileName);
    if (!file.is_open()) {
        // 可选：抛出异常或记录错误；此处静默失败
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // 跳过空行和注释
        if (line.empty() || line[0] == '#') continue;

        // 去除行首空白
        size_t start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) continue;
        line = line.substr(start);

        // 判断行类型（忽略大小写）
        char firstChar = static_cast<char>(std::tolower(static_cast<unsigned char>(line[0])));
        if (firstChar != 'v' && firstChar != 'f') continue;

        // 顶点: v x y z
        if (firstChar == 'v' && (line.size() == 1 || !std::isalpha(static_cast<unsigned char>(line[1])))) {
            std::istringstream iss(line.substr(1));
            double x, y, z;
            if (iss >> x >> y >> z) {
                points.push_back({ x, y, z });
            }
        }
        // 面: f v1 v2 v3 ...
        else if (firstChar == 'f') {
            std::istringstream iss(line.substr(1));
            std::vector<int> faceIndices;
            std::string token;

            while (iss >> token) {
                // 提取顶点索引（忽略纹理/法向）
                size_t slash = token.find('/');
                std::string idxStr = slash != std::string::npos ? token.substr(0, slash) : token;
                if (idxStr.empty()) continue;

                try {
                    int idx = std::stoi(idxStr);
                    if (idx < 0) {
                        // 负索引：-1 表示最后一个顶点
                        idx = static_cast<int>(points.size()) + idx + 1;
                    }
                    // OBJ 索引从 1 开始，我们保留 1-based
                    if (idx >= 1 && idx <= static_cast<int>(points.size())) {
                        faceIndices.push_back(idx);
                    }
                    else {
                        faceIndices.clear(); // 索引越界，丢弃整个面
                        break;
                    }
                }
                catch (...) {
                    faceIndices.clear();
                    break;
                }
            }

            // 仅处理三角形（严格 3 个顶点）
            if (faceIndices.size() == 3) {
                indices.insert(indices.end(), faceIndices.begin(), faceIndices.end());
                indices.push_back(0); // 用 0 分隔每个三角形
            }
            // 注：若需支持四边形等，可在此处添加三角剖分逻辑
        }
    }

    file.close();
}

void TriMesh::FixNormalsToOutside()
{
    if (points.empty() || indices.empty()) return;

    // === 1. 高效解析三角形 (1-based索引) ===
    std::vector<std::array<int, 3>> tris;
    for (size_t i = 0, start = 0; i <= indices.size(); ++i) {
        if (i == indices.size() || indices[i] == 0) {
            if (i - start == 3)
                tris.emplace_back(std::array{ indices[start], indices[start + 1], indices[start + 2] });
            start = i + 1;
        }
    }
    if (tris.empty()) return;

    // === 2. 构建边→三角形映射 (标准化边: min,max; 1-based索引存储避免0歧义) ===
    std::map<std::pair<int, int>, std::array<int, 2>> edgeMap;
    for (size_t i = 0; i < tris.size(); ++i) {
        const auto& t = tris[i];
        for (int j = 0; j < 3; ++j) {
            int a = t[j], b = t[(j + 1) % 3];
            if (a == b) continue; // 跳过退化边
            auto key = std::minmax(a, b);
            auto& entry = edgeMap[key];
            if (entry[0] == 0) entry[0] = static_cast<int>(i) + 1;
            else if (entry[1] == 0) entry[1] = static_cast<int>(i) + 1;
            // 非流形边（>2三角形共享）自动跳过：entry[1]被覆盖但后续BFS会因校验跳过
        }
    }

    // === 3. BFS统一局部法向 (处理所有连通分量) ===
    std::vector<bool> visited(tris.size(), false);
    auto needsFlip = [&](int t1_idx, int a, int b) -> bool {
        const auto& t1 = tris[t1_idx];
        for (int k = 0; k < 3; ++k) {
            if (t1[k] == a && t1[(k + 1) % 3] == b) return true;  // 顺序相同 → 法向相反
            if (t1[k] == b && t1[(k + 1) % 3] == a) return false; // 顺序相反 → 法向一致
        }
        return false; // 理论不可达（边映射已保证存在）
        };

    for (size_t seed = 0; seed < tris.size(); ++seed) {
        if (visited[seed]) continue;
        std::queue<int> q({ static_cast<int>(seed) });
        visited[seed] = true;

        while (!q.empty()) {
            int t0_idx = q.front(); q.pop();
            const auto& t0 = tris[t0_idx];

            for (int j = 0; j < 3; ++j) {
                int a = t0[j], b = t0[(j + 1) % 3];
                if (a == b) continue; // 防御性检查（构建时已过滤）
                auto it = edgeMap.find(std::minmax(a, b));
                if (it == edgeMap.end() || it->second[1] == 0) continue; // 非流形/无效边跳过

                // 获取相邻三角形索引（转0-based）
                int t1_idx = (it->second[0] - 1 == t0_idx) ? it->second[1] - 1 : it->second[0] - 1;
                if (visited[t1_idx]) continue;

                if (needsFlip(t1_idx, a, b))
                    std::swap(tris[t1_idx][1], tris[t1_idx][2]);

                visited[t1_idx] = true;
                q.push(t1_idx);
            }
        }
    }

    // === 4. 全局校正：有向体积决定内外朝向 ===
    double vol = 0.0;
    for (const auto& t : tris) {
        const auto& v0 = points[t[0] - 1], & v1 = points[t[1] - 1], & v2 = points[t[2] - 1];
        vol += v0.Dot(v1.Cross(v2)); // 6倍有向体积
    }
    if (vol < 0.0)
        for (auto& t : tris) std::swap(t[1], t[2]); // 全局翻转

    // === 5. 重建indices (1-based + 0分隔) ===
    indices.clear();
    for (const auto& t : tris)
        indices.insert(indices.end(), { t[0], t[1], t[2], 0 });
}

void TriMesh::CreateSphere(const Vec3d& center, double radius, int subdivisionLevel)
{
    // === 参数校验与初始化 ===
    if (radius <= std::numeric_limits<double>::epsilon()) {
        points.clear();
        indices.clear();
        return;
    }
    if (subdivisionLevel < 0) subdivisionLevel = 0;

    // 清空现有网格数据
    points.clear();
    indices.clear();

    // === 1. 构建初始八面体（单位球面，0-based索引）===
    std::vector<Vec3d> verts = {
        {1,0,0}, {- 1,0,0}, // X轴
        {0,1,0}, {0,-1,0}, // Y轴
        {0,0,1}, {0,0,-1}  // Z轴
    };

    // 初始8个三角形（经有向体积验证：法向严格向外）
    std::vector<std::array<int, 3>> faces = {
        {0,2,4}, {0,4,3}, {0,3,5}, {0,5,2}, // 包含顶点0的4个面
        {1,4,2}, {1,3,4}, {1,5,3}, {1,2,5}  // 包含顶点1的4个面
    };

    // === 2. 递归细分（每次将三角形分为4个）===
    for (int level = 0; level < subdivisionLevel; ++level) {
        std::map<std::pair<int, int>, int> edgeMap; // (minIdx, maxIdx) -> midpointIdx
        std::vector<std::array<int, 3>> newFaces;
        newFaces.reserve(faces.size() * 4);

        // 辅助：获取/创建边中点（归一化到单位球面）
        auto getMidpoint = [&](int i1, int i2) -> int {
            auto key = std::make_pair(std::min(i1, i2), std::max(i1, i2));
            if (auto it = edgeMap.find(key); it != edgeMap.end())
                return it->second;

            // 计算中点并归一化（若Vec3d无Normalization，替换为下方注释代码）
            Vec3d mid = (verts[i1] + verts[i2]) * 0.5;
            // double len = mid.Length(); if (len > 1e-12) mid /= len; // 手动归一化
            mid = mid.Normalization();
            int idx = static_cast<int>(verts.size());
            verts.push_back(mid);
            edgeMap[key] = idx;
            return idx;
            };

        // 细分每个三角形
        for (const auto& f : faces) {
            int a = f[0], b = f[1], c = f[2];
            int ab = getMidpoint(a, b);
            int bc = getMidpoint(b, c);
            int ca = getMidpoint(c, a);

            // 保持右手定则（法向向外）
            newFaces.push_back({ a, ab, ca });
            newFaces.push_back({ b, bc, ab });
            newFaces.push_back({ c, ca, bc });
            newFaces.push_back({ ab, bc, ca });
        }
        faces = std::move(newFaces);
    }

    // === 3. 应用缩放与平移（单位球 → 目标球）===
    points.reserve(verts.size());
    for (const Vec3d& v : verts) {
        points.push_back(center + v * radius);
    }

    // === 4. 构建TriMesh专用索引格式（1-based + 0分隔）===
    indices.reserve(faces.size() * 4);
    for (const auto& f : faces) {
        indices.push_back(f[0] + 1); // 转1-based
        indices.push_back(f[1] + 1);
        indices.push_back(f[2] + 1);
        indices.push_back(0);        // 面分隔符
    }
}

MESHMODELDLL double GetGlobalPrecision()
{
    return g_epsilon;
}

MESHMODELDLL void SetGlobalPrecision(double epsilon)
{
    g_epsilon = epsilon;
}
