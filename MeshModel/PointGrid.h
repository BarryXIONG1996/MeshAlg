#pragma once

#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <cstddef>

// 检测点维度
namespace detail {
    // 主模板: 默认情况(没有.z成员)
    template<typename T, typename = void>
    struct has_z_member : std::false_type {};

    // 偏特化模板: 当T由.z成员时,成功;否则,退化回主模板
    template<typename T>
    struct has_z_member<T, decltype(std::declval<T&>().z, void())> : std::true_type {};

    template<typename T>
    inline constexpr bool has_z_member_v = has_z_member<T>::value;

    template<typename T>
    struct point_dimension {
        static constexpr size_t value = has_z_member_v<T> ? 3 : 2;
    };

    template<typename T>
    inline constexpr size_t point_dimension_v = point_dimension<T>::value;
}

// 网格计算器模板
template<size_t Dim>
struct GridCalculator;

// 2D网格计算器
template<>
struct GridCalculator<2> {
    using Cell = std::array<int64_t, 2>;

    struct CellHash {
        size_t operator()(const Cell& c) const noexcept {
            size_t h = std::hash<int64_t>()(c[0]);
            h ^= std::hash<int64_t>()(c[1]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    static constexpr size_t neighbor_count = 9;
    static constexpr std::array<Cell, neighbor_count> neighbor_offsets = []() {
        std::array<Cell, neighbor_count> offsets{};
        size_t idx = 0;
        for (int64_t dx = -1; dx <= 1; ++dx)
            for (int64_t dy = -1; dy <= 1; ++dy)
                offsets[idx++] = { dx, dy };
        return offsets;
        }();

        template<typename Point>
        static Cell compute_cell(const Point& p, double cell_size_inv) noexcept {
            constexpr double MAX_SAFE = static_cast<double>(1LL << 53);
            auto safe_floor = [](double val) -> int64_t {
                if (!std::isfinite(val)) return 0;
                if (val >= MAX_SAFE) return std::numeric_limits<int64_t>::max() - 1;
                if (val <= -MAX_SAFE) return std::numeric_limits<int64_t>::min() + 1;
                return static_cast<int64_t>(std::floor(val));
                };
            return { safe_floor(p.x * cell_size_inv), safe_floor(p.y * cell_size_inv) };
        }

        template<typename Point1, typename Point2>
        static auto distance_sq(const Point1& a, const Point2& b) {
            auto dx = a.x - b.x;
            auto dy = a.y - b.y;
            return dx * dx + dy * dy;
        }

        template<typename Func>
        static void for_each_neighbor(const Cell& base_cell, Func&& func) {
            for (const auto& offset : neighbor_offsets) {
                Cell neighbor{ base_cell[0] + offset[0], base_cell[1] + offset[1] };
                func(neighbor);
            }
        }
};

// 3D网格计算器
template<>
struct GridCalculator<3> {
    using Cell = std::array<int64_t, 3>;

    struct CellHash {
        size_t operator()(const Cell& c) const noexcept {
            size_t h = std::hash<int64_t>()(c[0]);
            h ^= std::hash<int64_t>()(c[1]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int64_t>()(c[2]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    static constexpr size_t neighbor_count = 27;
    static constexpr std::array<Cell, neighbor_count> neighbor_offsets = []() {
        std::array<Cell, neighbor_count> offsets{};
        size_t idx = 0;
        for (int64_t dx = -1; dx <= 1; ++dx)
            for (int64_t dy = -1; dy <= 1; ++dy)
                for (int64_t dz = -1; dz <= 1; ++dz)
                    offsets[idx++] = { dx, dy, dz };
        return offsets;
        }();

        template<typename Point>
        static Cell compute_cell(const Point& p, double cell_size_inv) noexcept {
            constexpr double MAX_SAFE = static_cast<double>(1LL << 53);
            auto safe_floor = [](double val) -> int64_t {
                if (!std::isfinite(val)) return 0;
                if (val >= MAX_SAFE) return std::numeric_limits<int64_t>::max() - 1;
                if (val <= -MAX_SAFE) return std::numeric_limits<int64_t>::min() + 1;
                return static_cast<int64_t>(std::floor(val));
                };
            return {
                safe_floor(p.x * cell_size_inv),
                safe_floor(p.y * cell_size_inv),
                safe_floor(p.z * cell_size_inv)
            };
        }

        template<typename Point1, typename Point2>
        static auto distance_sq(const Point1& a, const Point2& b) {
            auto dx = a.x - b.x;
            auto dy = a.y - b.y;
            auto dz = a.z - b.z;
            return dx * dx + dy * dy + dz * dz;
        }

        template<typename Func>
        static void for_each_neighbor(const Cell& base_cell, Func&& func) {
            for (const auto& offset : neighbor_offsets) {
                Cell neighbor{ base_cell[0] + offset[0], base_cell[1] + offset[1], base_cell[2] + offset[2] };
                func(neighbor);
            }
        }
};

// 主模板
template<typename T, typename Point>
class PointGrid {
private:
    static constexpr size_t DIM = detail::point_dimension_v<Point>;
    static_assert(DIM == 2 || DIM == 3, "Point type must be 2D (x,y) or 3D (x,y,z)");

    using Calculator = GridCalculator<DIM>;
    using Cell = typename Calculator::Cell;
    using CellHash = typename Calculator::CellHash;
    using StorageVector = std::vector<std::pair<Point, T>>;
    using GridMap = std::unordered_map<Cell, StorageVector, CellHash>;

    GridMap grid_;
    double cell_size_;
    double cell_size_inv_;
    double eps_sq_;
    size_t total_size_ = 0;

    // 前向声明迭代器类
    class const_iterator;
    class iterator;

    // 核心插入逻辑
    template<typename Key, typename Val>
    std::pair<iterator, bool> insert_impl(Key&& key, Val&& val) {
        if (auto it = find(key); it != end())
            return { it, false };

        Cell cell = Calculator::compute_cell(key, cell_size_inv_);
        auto& vec = grid_[cell];
        auto vec_it = vec.emplace(
            vec.end(),
            std::forward<Key>(key),
            std::forward<Val>(val)
        );
        ++total_size_;
        return { iterator(grid_.find(cell), vec_it, this), true };
    }

    void initialize(double epsilon) {
        if (epsilon <= 0)
            throw std::invalid_argument("epsilon must be positive");
        constexpr double MIN_EPS = 1e-12;
        cell_size_ = std::max(epsilon, MIN_EPS);
        cell_size_inv_ = 1.0 / cell_size_;
        eps_sq_ = cell_size_ * cell_size_;
    }

public:
    // 在public部分重新声明迭代器类型
    class const_iterator;
    class iterator;

    using key_type = Point;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    // ============ 迭代器实现 ============
    class const_iterator {
    private:
        using GridIter = typename GridMap::const_iterator;
        using VecIter = typename StorageVector::const_iterator;

        GridIter grid_it_;
        GridIter grid_end_;
        VecIter vec_it_;
        const PointGrid* parent_;

        void advance() noexcept {
            while (grid_it_ != grid_end_ && vec_it_ == grid_it_->second.end()) {
                ++grid_it_;
                if (grid_it_ != grid_end_)
                    vec_it_ = grid_it_->second.begin();
            }
        }

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = PointGrid::value_type;
        using difference_type = PointGrid::difference_type;
        using pointer = const value_type*;
        using reference = const value_type&;

        const_iterator() noexcept : parent_(nullptr) {}

        // 默认构造的迭代器指向end
        const_iterator(GridIter g, VecIter v, const PointGrid* p) noexcept
            : grid_it_(g), grid_end_(p ? p->grid_.end() : GridIter{}),
            vec_it_(v), parent_(p) {
            if (parent_ && grid_it_ != grid_end_) advance();
        }

        reference operator*() const noexcept {
            if (parent_ == nullptr)
                throw std::runtime_error("Dereferencing end iterator");
            return *vec_it_;
        }

        pointer operator->() const noexcept {
            if (parent_ == nullptr)
                throw std::runtime_error("Dereferencing end iterator");
            return &*vec_it_;
        }

        const_iterator& operator++() {
            if (parent_ && grid_it_ != grid_end_) {
                ++vec_it_;
                advance();
            }
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const const_iterator& o) const {
            // 安全性检查：如果parent_不同，则迭代器必然不相等
            if (parent_ != o.parent_) return false;
            // 如果都是nullptr，都是默认构造的end迭代器
            if (parent_ == nullptr) return true;
            // 否则比较内部迭代器
            if (grid_it_ != o.grid_it_) return false;
            if (grid_it_ == grid_end_) return true;
            return vec_it_ == o.vec_it_;
        }

        bool operator!=(const const_iterator& o) const { return !(*this == o); }
    };

    class iterator {
    private:
        using GridIter = typename GridMap::iterator;
        using VecIter = typename StorageVector::iterator;

        GridIter grid_it_;
        GridIter grid_end_;
        VecIter vec_it_;
        PointGrid* parent_;

        void advance() noexcept {
            while (grid_it_ != grid_end_ && vec_it_ == grid_it_->second.end()) {
                ++grid_it_;
                if (grid_it_ != grid_end_)
                    vec_it_ = grid_it_->second.begin();
            }
        }

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = PointGrid::value_type;
        using difference_type = PointGrid::difference_type;
        using pointer = value_type*;
        using reference = value_type&;

        iterator() noexcept : parent_(nullptr) {}

        iterator(GridIter g, VecIter v, PointGrid* p) noexcept
            : grid_it_(g), grid_end_(p ? p->grid_.end() : GridIter{}),
            vec_it_(v), parent_(p) {
            if (parent_ && grid_it_ != grid_end_) advance();
        }

        reference operator*() const noexcept {
            if (parent_ == nullptr)
                throw std::runtime_error("Dereferencing end iterator");
            return *vec_it_;
        }

        pointer operator->() const noexcept {
            if (parent_ == nullptr)
                throw std::runtime_error("Dereferencing end iterator");
            return &*vec_it_;
        }

        iterator& operator++() {
            if (parent_ && grid_it_ != grid_end_) {
                ++vec_it_;
                advance();
            }
            return *this;
        }

        iterator operator++(int) {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const iterator& o) const {
            // 安全性检查：如果parent_不同，则迭代器必然不相等
            if (parent_ != o.parent_) return false;
            // 如果都是nullptr，都是默认构造的end迭代器
            if (parent_ == nullptr) return true;
            // 否则比较内部迭代器
            if (grid_it_ != o.grid_it_) return false;
            if (grid_it_ == grid_end_) return true;
            return vec_it_ == o.vec_it_;
        }

        bool operator!=(const iterator& o) const { return !(*this == o); }

        operator const_iterator() const {
            return const_iterator(grid_it_, vec_it_, parent_);
        }
    };

    // ============ 构造/析构 ============
    explicit PointGrid(double eps = 1e-6) { initialize(eps); }
    PointGrid(const PointGrid&) = default;
    PointGrid(PointGrid&&) noexcept = default;
    PointGrid& operator=(const PointGrid&) = default;
    PointGrid& operator=(PointGrid&&) noexcept = default;
    ~PointGrid() = default;

    // ============ 容量 ============
    bool empty() const noexcept { return total_size_ == 0; }
    size_type size() const noexcept { return total_size_; }
    void clear() noexcept {
        grid_.clear();
        total_size_ = 0;
    }
    void reserve(size_type n) { grid_.reserve(n); }
    size_type bucket_count() const noexcept { return grid_.bucket_count(); }

    // ============ 查找 ============
    iterator find(const key_type& key) {
        Cell base = Calculator::compute_cell(key, cell_size_inv_);

        for (const auto& offset : Calculator::neighbor_offsets) {
            Cell nb = base;
            for (size_t i = 0; i < DIM; ++i) nb[i] += offset[i];

            auto cell_it = grid_.find(nb);
            if (cell_it != grid_.end()) {
                for (auto vec_it = cell_it->second.begin(); vec_it != cell_it->second.end(); ++vec_it) {
                    if (Calculator::distance_sq(vec_it->first, key) < eps_sq_) {
                        return iterator(cell_it, vec_it, this); // 返回找到的第一个元素
                    }
                }
            }
        }
        return end();
    }

    const_iterator find(const key_type& key) const {
        Cell base = Calculator::compute_cell(key, cell_size_inv_);

        for (const auto& offset : Calculator::neighbor_offsets) {
            Cell nb = base;
            for (size_t i = 0; i < DIM; ++i) nb[i] += offset[i];

            auto cell_it = grid_.find(nb);
            if (cell_it != grid_.end()) {
                for (auto vec_it = cell_it->second.begin(); vec_it != cell_it->second.end(); ++vec_it) {
                    if (Calculator::distance_sq(vec_it->first, key) < eps_sq_) {
                        return const_iterator(cell_it, vec_it, this);
                    }
                }
            }
        }
        return cend();
    }

    bool contains(const key_type& key) const { return find(key) != cend(); }
    size_type count(const key_type& key) const { return contains(key) ? 1 : 0; }

    // ============ 访问 ============
    const mapped_type& at(const key_type& key) const {
        if (auto it = find(key); it != cend())
            return it->second;
        throw std::out_of_range("PointGrid::at: key not found");
    }

    mapped_type& at(const key_type& key) {
        if (auto it = find(key); it != end())
            return it->second;
        throw std::out_of_range("PointGrid::at: key not found");
    }

    mapped_type& operator[](const key_type& key) {
        auto [it, inserted] = emplace(key, mapped_type{});
        return it->second;
    }

    const mapped_type& operator[](const key_type& key) const { return at(key); }

    // ============ 插入 ============
    std::pair<iterator, bool> insert(const value_type& val) {
        return insert_impl(val.first, val.second);
    }

    std::pair<iterator, bool> insert(value_type&& val) {
        return insert_impl(val.first, std::move(val.second));
    }

    template<typename P,
        std::enable_if_t<std::is_constructible_v<value_type, P&&>, int> = 0>
    std::pair<iterator, bool> insert(P&& val) {
        return insert_impl(
            static_cast<const key_type&>(std::forward<P>(val).first),
            std::forward<P>(val).second
        );
    }

    void insert(std::initializer_list<value_type> ilist) {
        for (const auto& v : ilist) insert(v);
    }

    template<typename... Args>
    std::pair<iterator, bool> emplace(const key_type& key, Args&&... args) {
        return insert_impl(key, mapped_type(std::forward<Args>(args)...));
    }

    template<typename... Args>
    std::pair<iterator, bool> emplace(key_type&& key, Args&&... args) {
        return insert_impl(std::move(key), mapped_type(std::forward<Args>(args)...));
    }

    template<typename M>
    std::pair<iterator, bool> insert_or_assign(const key_type& key, M&& obj) {
        if (auto it = find(key); it != end()) {
            it->second = std::forward<M>(obj);
            return { it, false };
        }
        return emplace(key, std::forward<M>(obj));
    }

    // ============ 删除 ============
    size_type erase(const key_type& key) {
        Cell base = Calculator::compute_cell(key, cell_size_inv_);

        for (const auto& offset : Calculator::neighbor_offsets) {
            Cell nb = base;
            for (size_t i = 0; i < DIM; ++i) nb[i] += offset[i];

            auto cell_it = grid_.find(nb);
            if (cell_it != grid_.end()) {
                auto& vec = cell_it->second;
                for (auto it = vec.begin(); it != vec.end(); ++it) {
                    if (Calculator::distance_sq(it->first, key) < eps_sq_) {
                        vec.erase(it);
                        if (vec.empty()) grid_.erase(cell_it);
                        --total_size_;
                        return 1;
                    }
                }
            }
        }
        return 0;
    }

     iterator erase(iterator pos) {
        if (pos == end()) return end();
    
        auto cell_it = pos.grid_it_;
        auto& vec = cell_it->second;
    
        // erase 返回下一个有效 vec iterator
        auto next_vec_it = vec.erase(pos.vec_it_);
        --total_size_;
    
        if (vec.empty()) {
            // 整个 cell 为空，从 grid_ 中移除
            auto next_cell_it = grid_.erase(cell_it);
            // 构造指向下一个非空 cell 的 begin（如果存在）
            if (next_cell_it != grid_.end()) {
                return iterator(next_cell_it, next_cell_it->second.begin(), this);
            } else {
                return end();
            }
        } else {
            // 同一个 cell 内还有元素
            return iterator(cell_it, next_vec_it, this);
        }
    }

    iterator erase(iterator first, iterator last) {
        if (first == last) return first;
        std::vector<key_type> keys;
        for (auto it = first; it != last; ++it)
            keys.push_back(it->first);
        for (const auto& k : keys)
            erase(k);
        return end();
    }

    // ============ 迭代器 ============
    iterator begin() {
        auto it = grid_.begin();
        while (it != grid_.end() && it->second.empty()) ++it;
        if (it != grid_.end())
            return iterator(it, it->second.begin(), this);
        return end();
    }

    iterator end() {
        return iterator(grid_.end(), typename StorageVector::iterator{}, this);
    }

    const_iterator begin() const {
        auto it = grid_.begin();
        while (it != grid_.end() && it->second.empty()) ++it;
        if (it != grid_.end())
            return const_iterator(it, it->second.begin(), this);
        return cend();
    }

    const_iterator end() const {
        return const_iterator(grid_.end(), typename StorageVector::const_iterator{}, this);
    }

    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

    // ============ 信息查询 ============
    double epsilon() const noexcept { return cell_size_; }
    size_type bucket_size(size_type n) const {
        if (n >= bucket_count())
            throw std::out_of_range("PointGrid::bucket_size: index out of range");
        return grid_.bucket_size(n);
    }
    double load_factor() const noexcept { return grid_.load_factor(); }
    double max_load_factor() const noexcept { return grid_.max_load_factor(); }
    void max_load_factor(float ml) { grid_.max_load_factor(ml); }
    void rehash(size_type n) { grid_.rehash(n); }

    // ============ 交换 ============
    void swap(PointGrid& other) {
        const double max_eps = std::max(cell_size_, other.cell_size_);
        if (std::abs(cell_size_ - other.cell_size_) > max_eps * 1e-6)
            throw std::invalid_argument("PointGrid::swap: epsilon mismatch");
        std::swap(grid_, other.grid_);
        std::swap(cell_size_, other.cell_size_);
        std::swap(cell_size_inv_, other.cell_size_inv_);
        std::swap(eps_sq_, other.eps_sq_);
        std::swap(total_size_, other.total_size_);
    }
};

// 非成员函数
template<typename T, typename Point>
void swap(PointGrid<T, Point>& a, PointGrid<T, Point>& b) {
    a.swap(b);
}

template<typename T, typename Point>
bool operator==(const PointGrid<T, Point>& lhs, const PointGrid<T, Point>& rhs) {
    if (lhs.size() != rhs.size()) return false;
    const double max_eps = std::max(lhs.epsilon(), rhs.epsilon());
    if (std::abs(lhs.epsilon() - rhs.epsilon()) > max_eps * 1e-6) return false;
    for (const auto& elem : lhs) {
        auto it = rhs.find(elem.first);
        if (it == rhs.cend() || !(it->second == elem.second))
            return false;
    }
    return true;
}

template<typename T, typename Point>
bool operator!=(const PointGrid<T, Point>& lhs, const PointGrid<T, Point>& rhs) {
    return !(lhs == rhs);
}
