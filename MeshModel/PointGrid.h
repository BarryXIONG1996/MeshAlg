#pragma once
#ifndef POINT_GRID_H
#define POINT_GRID_H

#include <unordered_map>
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <cassert>
#include <cstddef>
#include "Geometry.h" // 包含 Vec3d 定义

#ifdef MESHMODELDLL_EXPORTS
#define MESHMODELDLL __declspec(dllexport)
#else
#define MESHMODELDLL __declspec(dllimport)
#endif

/**
 * @brief 网格哈希点去重器（std::map 接口兼容版｜工业级安全）
 * @tparam T 关联值类型（Vertex*, size_t, shared_ptr 等）
 *
 * @核心设计:
 * 1. 完全兼容 std::map<Vec3d, T, Cmp> 的公共接口（insert/find/count/operator[]/at）
 * 2. 严格去重语义：‖P-Q‖ < epsilon -> 视为同一点（27邻域+精确距离验证）
 * 3. 安全防护：所有返回引用均为 const，禁止通过迭代器/[] 修改内部值
 * 4. 明确行为：operator[] 不存在时抛异常（拒绝隐式插入），符合去重逻辑
 *
 * @与 std::map 关键差异说明（文档强制要求）:
 * | 操作 | std::map 行为 | PointGrid 行为 | 原因 |
 * |------|----------------|-----------------|------|
 * | operator[] | 不存在时插入 T() | 不存在时抛 std::out_of_range | 避免破坏去重逻辑（隐式插入默认值会导致后续近似点被错误合并） |
 * | insert | 总是插入新元素 | 近似点存在时返回已存在元素+false | 实现几何去重核心语义 |
 * | 迭代器解引用 | pair<const Key, T>& | pair<const Vec3d, T>& | Key 语义为 const（坐标不应被修改） |
 * | 修改值 | 允许 p->second = ... | 禁止（返回 const 引用） | 防止破坏去重拓扑一致性 |
 */
template<typename T>
class MESHMODELDLL PointGrid
{
private:
    // ====== 类型别名（严格对标 std::map） ======
    using key_type = Vec3d;
    using mapped_type = T;
    using value_type = std::pair<const Vec3d, T>; // 注意：first 为 const，符合 STL 习惯
    using size_type = size_t;

    // ====== 私有辅助结构（MSVC 兼容：置于公有前） ======
    struct Cell {
        int x, y, z;
        bool operator==(const Cell& o) const noexcept {
            return x == o.x && y == o.y && z == o.z;
        }
    };

    struct CellHash {
        size_t operator()(const Cell& c) const noexcept {
            size_t h1 = std::hash<int>()(c.x);
            size_t h2 = std::hash<int>()(c.y);
            size_t h3 = std::hash<int>()(c.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h1 >> 3);
        }
    };

    // ====== 存储结构（关键：存储 pair<const Vec3d, T>） ======
    // 为满足 value_type = pair<const Vec3d, T>，内部使用 pair<Vec3d, T> 存储，
    // 但迭代器/接口层强制视为 const（通过类型转换安全实现）
    using StoragePair = std::pair<Vec3d, T>;
    std::unordered_map<Cell, std::vector<StoragePair>, CellHash> grid_;

    const double cell_size_;
    const double cell_size_inv_;
    const double eps_sq_;

    // ====== 私有工具函数 ======
    Cell ComputeCell(const Vec3d& p) const noexcept {
        return {
            static_cast<int>(std::floor(p.x * cell_size_inv_)),
            static_cast<int>(std::floor(p.y * cell_size_inv_)),
            static_cast<int>(std::floor(p.z * cell_size_inv_))
        };
    }

    // 安全转换：StoragePair* -> value_type*（仅用于迭代器解引用）
    static const value_type* ToValueTypePtr(const StoragePair* p) noexcept {
        // reinterpret_cast 安全：仅改变 first 的 const 限定符，内存布局完全相同
        return reinterpret_cast<const value_type*>(p);
    }

public:
    // ==================== 迭代器（严格对标 std::map::const_iterator） ====================
    class const_iterator {
        friend class PointGrid;
        using GridIt = typename decltype(grid_)::const_iterator;
        using CellVecIt = typename std::vector<StoragePair>::const_iterator;

        GridIt outer_;
        CellVecIt inner_;
        const PointGrid* parent_;

        const_iterator(GridIt o, CellVecIt i, const PointGrid* p)
            : outer_(o), inner_(i), parent_(p) {}

    public:
        using value_type = PointGrid::value_type; // pair<const Vec3d, T>
        using reference = const value_type&;
        using pointer = const value_type*;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        reference operator*() const noexcept {
            return *PointGrid::ToValueTypePtr(&*inner_);
        }
        pointer operator->() const noexcept {
            return PointGrid::ToValueTypePtr(&*inner_);
        }

        const_iterator& operator++() {
            if (outer_ == parent_->grid_.end()) return *this;
            ++inner_;
            while (outer_ != parent_->grid_.end() && inner_ == outer_->second.end()) {
                ++outer_;
                if (outer_ != parent_->grid_.end()) inner_ = outer_->second.begin();
            }
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const const_iterator& o) const noexcept {
            if (parent_ != o.parent_) return false;
            bool end1 = (outer_ == parent_->grid_.end());
            bool end2 = (o.outer_ == o.parent_->grid_.end());
            if (end1 || end2) return end1 == end2;
            return outer_ == o.outer_ && inner_ == o.inner_;
        }

        bool operator!=(const const_iterator& o) const noexcept {
            return !(*this == o);
        }
    };

    using iterator = const_iterator; // 仅提供 const 迭代器（禁止修改）

    // ==================== 构造/析构 ====================
    explicit PointGrid(double epsilon)
        : cell_size_(epsilon > 1e-12 ? epsilon : 1e-6)
        , cell_size_inv_(1.0 / cell_size_)
        , eps_sq_(cell_size_* cell_size_)
    {
        assert(epsilon > 0 && "epsilon must be positive");
    }

    ~PointGrid() = default;
    PointGrid(const PointGrid&) = delete;
    PointGrid& operator=(const PointGrid&) = delete;

    // ==================== std::map 核心接口（完全兼容） ====================

    /**
     * @brief 插入元素（去重核心）
     * @param value std::pair<Vec3d, T>（注意：first 非 const，但存储后视为 const）
     * @return pair<iterator, bool>
     *   - bool=true: 新插入，iterator 指向新元素
     *   - bool=false: 存在近似点，iterator 指向已存在元素（值未被覆盖！）
     * @note 与 std::map::insert 行为差异：
     *       std::map 总是插入新元素；PointGrid 在近似点存在时拒绝插入并返回已有元素
     */
    std::pair<iterator, bool> insert(const value_type& value) {
        return insert_impl(value.first, value.second);
    }

    std::pair<iterator, bool> insert(value_type&& value) {
        return insert_impl(std::move(value.first), std::move(value.second));
    }

    /**
     * @brief 查询元素（对标 std::map::find）
     * @return 指向元素的迭代器；不存在则返回 end()
     * @note 与 std::map::find 行为一致（仅查询，不插入）
     */
    iterator find(const key_type& key) const noexcept {
        Cell base = ComputeCell(key);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Cell nb{ base.x + dx, base.y + dy, base.z + dz };
                    auto cell_it = grid_.find(nb);
                    if (cell_it == grid_.end()) continue;

                    const auto& cell_vec = cell_it->second;
                    for (auto vec_it = cell_vec.begin(); vec_it != cell_vec.end(); ++vec_it) {
                        if ((vec_it->first - key).LengthSquared() < eps_sq_) {
                            return iterator(cell_it, vec_it, this);
                        }
                    }
                }
            }
        }
        return end();
    }

    /**
     * @brief 计数（对标 std::map::count）
     * @return 0 或 1（去重保证唯一点）
     */
    size_type count(const key_type& key) const noexcept {
        return find(key) != end() ? 1 : 0;
    }

    /**
     * @brief 安全访问（对标 std::map::at）
     * @return 元素的 const 引用；不存在则抛 std::out_of_range
     */
    const mapped_type& at(const key_type& key) const {
        auto it = find(key);
        if (it == end()) {
            throw std::out_of_range("PointGrid::at: key not found (within epsilon)");
        }
        return it->second;
    }

    /**
     * @brief 下标访问（对标 std::map::operator[]，但行为关键差异！）
     * @return 元素的 const 引用
     * @throw std::out_of_range 若点不存在（拒绝隐式插入！）
     * @warning 与 std::map 的核心差异：
     *          std::map 在 key 不存在时插入 T()；PointGrid 严格抛异常
     *          原因：隐式插入默认值会破坏去重逻辑（后续近似点将错误合并到默认值）
     *          使用前务必确认点已通过 insert() 插入！
     */
    const mapped_type& operator[](const key_type& key) const {
        return at(key);
    }

    // 非 const 对象调用也返回 const 引用（禁止修改）
    const mapped_type& operator[](const key_type& key) {
        return const_cast<const PointGrid*>(this)->at(key);
    }

    // ==================== 迭代器接口（完全对标 std::map） ====================
    iterator begin() const noexcept {
        if (grid_.empty()) return end();
        auto it = grid_.begin();
        return iterator(it, it->second.begin(), this);
    }

    iterator end() const noexcept {
        return iterator(grid_.end(), {}, this);
    }

    iterator cbegin() const noexcept { return begin(); }
    iterator cend() const noexcept { return end(); }

    // ==================== 辅助接口 ====================
    void clear() noexcept { grid_.clear(); }          // 小写：对标 STL
    bool empty() const noexcept { return grid_.empty(); }
    size_type size() const noexcept {                 // 对标 STL
        size_type cnt = 0;
        for (const auto& cell : grid_) cnt += cell.second.size();
        return cnt;
    }
    double epsilon() const noexcept { return cell_size_; } // 小写：符合 STL 命名习惯

private:
    // insert 的实现核心（避免代码重复）
    template<typename Key, typename Val>
    std::pair<iterator, bool> insert_impl(Key&& key, Val&& val) {
        Cell base = ComputeCell(key);

        // 检查 27 邻域（含自身）
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Cell nb{ base.x + dx, base.y + dy, base.z + dz };
                    auto cell_it = grid_.find(nb);
                    if (cell_it == grid_.end()) continue;

                    auto& cell_vec = cell_it->second;
                    for (auto vec_it = cell_vec.begin(); vec_it != cell_vec.end(); ++vec_it) {
                        if ((vec_it->first - key).LengthSquared() < eps_sq_) {
                            // 找到近似点：返回已有元素迭代器 + false（拒绝插入）
                            return { iterator(cell_it, vec_it, this), false };
                        }
                    }
                }
            }
        }

        // 无近似点：插入新元素
        auto cell_it = grid_.find(base);
        if (cell_it == grid_.end()) {
            cell_it = grid_.emplace(base, std::vector<StoragePair>()).first;
        }
        auto& cell_vec = cell_it->second;
        size_t idx = cell_vec.size();
        cell_vec.emplace_back(std::forward<Key>(key), std::forward<Val>(val));
        return { iterator(cell_it, cell_vec.begin() + idx, this), true };
    }
};

// ==================== 显式实例化声明（按需启用） ====================
#ifdef MESHMODELDLL_EXPORTS
extern template class MESHMODELDLL PointGrid<Vertex*>;
extern template class MESHMODELDLL PointGrid<size_t>;
extern template class MESHMODELDLL PointGrid<int>;
#endif

#endif // POINT_GRID_H