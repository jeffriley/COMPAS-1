#pragma once


#include <cstddef>   // std::ptrdiff_t
#include <iterator>  // std::iterator_traits
#include <memory>    // std::pointer_traits, std::to_address


// Generic, reusable C++ alternative to boost::indirect_iterator.
// The functionality is equivalent - just doesn't require boost.
// Dereferences through whatever pointer-like element the base iterator yields, as const.
template <typename BaseIter>
class ConstIndirectIter {
    BaseIter m_It;
    using Ptr  = typename std::iterator_traits<BaseIter>::value_type;   // e.g. std::unique_ptr<T>
    using Elem = typename std::pointer_traits<Ptr>::element_type;       // T  (deduced, never named)
public:
    using iterator_category = typename std::iterator_traits<BaseIter>::iterator_category;
    using value_type        = Elem;
    using reference         = const Elem&;
    using pointer           = const Elem*;
    using difference_type   = std::ptrdiff_t;

    explicit ConstIndirectIter(BaseIter p_It) : m_It(p_It) {}

    const Elem& operator*()  const { return *(*m_It); }                 // const ptr& -> Elem& -> const Elem&
    const Elem* operator->() const { return std::to_address(*m_It); }   // works for unique_ptr/shared_ptr/raw
    ConstIndirectIter& operator++() { ++m_It; return *this; }
    ConstIndirectIter& operator--() { --m_It; return *this; }
    bool operator==(const ConstIndirectIter& o) const { return m_It == o.m_It; }
    bool operator!=(const ConstIndirectIter& o) const { return m_It != o.m_It; }
};
