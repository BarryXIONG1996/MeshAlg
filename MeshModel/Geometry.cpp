#include "Geometry.h"

#include <fstream>
#include <sstream>
#include <cctype>
#include <array>
#include <queue>
#include <cmath>
#include <algorithm>
#include <numeric>

const double g_epsilon = 1e-9;

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

bool Vec3dCmp::operator()(const Vec3d& a, const Vec3d& b) const {
    if (fabs(a.x - b.x) > g_epsilon) return a.x < b.x;
    if (fabs(a.y - b.y) > g_epsilon) return a.y < b.y;
    return a.z < b.z - g_epsilon;
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