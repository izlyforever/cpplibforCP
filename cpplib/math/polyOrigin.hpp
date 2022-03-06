#pragma once
#include "mod.hpp"
#include "poly.hpp"

template<typename T>
class PolyBaseOrigin : public PolyBase<T> {
 public:
  using PolyBase<T>::PolyBase;
  PolyBaseOrigin (const PolyBase<T>& x) : PolyBase<T>(x) {}
  PolyBaseOrigin (PolyBase<T>&& x) : PolyBase<T>(std::forward<PolyBase<T>>(x)) {}
 protected:
  PolyBaseOrigin mul(const PolyBaseOrigin& rhs) const {
    std::vector<T> ans(this->size() + rhs.size() - 1);
    for (int i = 0, sn = (int)this->size(); i < sn; ++i) {
      for (int j = 0, rsn = rhs.size(); j < rsn; ++j) {
        ans[i + j] += (*this)[i] * rhs[j];
      }
    }
    return PolyBaseOrigin(ans);
  }
};

// Origin Poly used for testing
constexpr int ORGM = 1e9 + 7;
using PolyOrigin = Poly<PolyBaseOrigin<MInt<ORGM>>, MInt<ORGM>>;
using PolyOriginDynamic = Poly<PolyBaseOrigin<ModInt>, ModInt>;
