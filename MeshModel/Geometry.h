#pragma once

#ifdef MESHMODEL_EXPORTS
#define MESHMODELDLL __declspec(dllexport)
#else
#define MESHMODELDLL __declspec(dllimport)
#endif

#include <vector>
#include <string>

MESHMODELDLL double GetGlobalPrecision();
MESHMODELDLL void SetGlobalPrecision(double epsilon);

struct Vec3d;
/**
 * @brief 4x4 双精度变换矩阵（列主序存储）
 * @note  所有函数定义位于 Mat4d.cpp
 *        内存布局: m[col*4 + row] (OpenGL/Vulkan 兼容)
 */
struct MESHMODELDLL Mat4d
{
    double m[16]; // 列主序: m[0]=m00, m[4]=m01, m[1]=m10...

    // ========== 构造函数声明 ==========
    Mat4d();
    explicit Mat4d(const double* data);
    Mat4d(double m00, double m01, double m02, double m03,
        double m10, double m11, double m12, double m13,
        double m20, double m21, double m22, double m23,
        double m30, double m31, double m32, double m33);

    // ========== 静态工厂函数声明 ==========
    static Mat4d Identity();
    static Mat4d Translation(const Vec3d& t);
    static Mat4d Scaling(const Vec3d& s);
    static Mat4d Rotation(const Vec3d& axis, double angle);
    static Mat4d RotationX(double angle);
    static Mat4d RotationY(double angle);
    static Mat4d RotationZ(double angle);

    // ========== 成员函数声明 ==========
    void SetIdentity();
    void Transpose();
    Mat4d operator*(const Mat4d& other) const;
    Vec3d TransformPoint(const Vec3d& v) const;
    Vec3d TransformVector(const Vec3d& v) const;
    double& operator()(int row, int col);
    double operator()(int row, int col) const;
};

struct MESHMODELDLL Vec2d
{
    double x, y;
    Vec2d operator+(const Vec2d& v) const;
    Vec2d operator-(const Vec2d& v) const;
    Vec2d operator*(double s) const;
    Vec2d operator/(double s) const;
    double Dot(const Vec2d& v) const;
    double Cross(const Vec2d& v) const;
    double Length() const;
    double LengthSq() const;
    Vec2d Normalization() const;
    bool Parallel(const Vec2d& v) const;
    bool Equal(const Vec2d& v) const;
};

struct MESHMODELDLL Vec3d
{
    double x, y, z;
    Vec3d operator-(const Vec3d& v) const;
    Vec3d operator+(const Vec3d& v) const;
    Vec3d operator*(double s) const;
    Vec3d operator/(double s) const;
    double Dot(const Vec3d& v) const;
    Vec3d Cross(const Vec3d& v) const;
    double Length() const;
    double LengthSq() const;
    Vec3d Normalization() const;
    bool Parallel(const Vec3d& v) const;
    bool Equal(const Vec3d& v) const;
};

struct MESHMODELDLL BndBox3d
{
    BndBox3d();
    Vec3d lowerBnd, upperBnd;
    BndBox3d Intersect(BndBox3d const& bnd) const;
    bool IsOut(BndBox3d bnd) const;
    void Add(const Vec3d& p);
    Vec3d Center() const;
};

struct MESHMODELDLL TriMesh
{
    std::vector<int> indices;
    std::vector<Vec3d> points;

    void Transform(Mat4d const& mat4d);
    void BuildFromOBJ(const std::string& fileName);
    void FixNormalsToOutside(); // 封闭实体法向修正
    void CreateSphere(const Vec3d& center, double radius, int subdivisionLevel); // nFaces = 8 * 4^n
};