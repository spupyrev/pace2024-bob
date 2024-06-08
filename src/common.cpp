#include "common.h"

#include <random>
#include <limits>
#include <cassert>

using namespace std;

std::vector<string> SplitNotNull(const std::string& ss, const std::string& c) {
  std::string s = ss + c;
  std::vector<std::string> result;
  std::string cur = "";

  for (size_t i = 0; i < s.length(); i++) {
    if (c.find(s[i]) != std::string::npos) {
      if (cur.length() > 0)
        result.push_back(cur);
      cur = "";
    } else {
      cur += s[i];
    }
  }
  return result;
}

std::vector<int> SplitNotNullInt(const string& ss, const string& c) {
  auto tmp = SplitNotNull(ss, c);
  vector<int> res;
  for (auto v : tmp) {
    res.push_back(to_int(v));
  }
  return res;
}

double Sum(const vector<double>& vec) {
  double sum = 0;
  for (double v : vec) {
    sum += v;
  }
  return sum;
}

double Sum(const vector<uint32_t>& vec) {
  double sum = 0;
  for (uint32_t v : vec) {
    sum += v;
  }
  return sum;
}

double Percentile(const vector<double>& vec, int p) {
  if (vec.empty())
    return 0;

  size_t n = vec.size();
  size_t pos = std::min(p * n / 100, n - 1);
  return vec[pos];
}

int Compare(double numberA, double numberB) {
  double diff = numberA - numberB;

  if (diff <= -EPS)
    return -1;
  if (diff >= EPS)
    return 1;
  return 0;
}

bool Equal(double a, double b) {
  return Abs(a - b) <= EPS;
}

bool Greater(double numberA, double numberB) {
  return Compare(numberA, numberB) > 0;
}

bool GreaterOrEqual(double numberA, double numberB) {
  return Compare(numberA, numberB) >= 0;
}

bool LessOrEqual(double numberA, double numberB) {
  return Compare(numberA, numberB) <= 0;
}

bool Less(double numberA, double numberB) {
  return Compare(numberA, numberB) < 0;
}

//std::mt19937 RNG;

// size_t Rand::setSeed() {
//   return setSeed(static_cast<size_t>(time(0)));
// }

// size_t Rand::setSeed(size_t seed) {
//   RNG.seed(seed);
//   return seed;
// }

// double Rand::nextDouble() {
//   return nextDouble(RNG);
// }

double Rand::nextDouble(std::mt19937& rng) {
  return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
}

// bool Rand::check(double probability) {
//   return nextDouble() <= probability;
// }

bool Rand::check(double probability, std::mt19937& rng) {
  return nextDouble(rng) <= probability;
}

// int Rand::next() {
//   return RNG();
// }

// int Rand::next(int bound) {
//   return Abs(RNG()) % bound;
// }

int Rand::next(int bound, std::mt19937& rng) {
  return Abs(rng()) % bound;
}

// int Rand::next(int lower, int upper) {
//   assert(0 <= lower && lower < upper);
//   return lower + Abs(RNG()) % (upper - lower);
// }

// std::vector<uint32_t> Rand::permutation(uint32_t n) {
//   return permutation(n, RNG);
// }

std::vector<uint32_t> Rand::permutation(uint32_t n, std::mt19937& rng) {
  std::vector<uint32_t> p(n);
  for (uint32_t i = 0; i < n; i++) {
    p[i] = i;
  }
  Rand::shuffle(p.begin(), p.end(), rng);
  return p;
}

std::vector<uint32_t> identity(uint32_t n) {
  std::vector<uint32_t> p(n);
  for (uint32_t i = 0; i < n; i++) {
    p[i] = i;
  }
  return p;
}

std::vector<uint32_t> inverse(const std::vector<uint32_t>& p) {
  const uint32_t NOT_SET = uint32_t(-1);
  std::vector<uint32_t> r(p.size(), NOT_SET);
  for (uint32_t i = 0; i < p.size(); i++) {
    assert(r[p[i]] == NOT_SET && "incorrect permutation for inversion");
    r[p[i]] = i;
  }
  return r;
}

std::vector<uint32_t> reverse(const std::vector<uint32_t>& p) {
  std::vector<uint32_t> r(p.begin(), p.end());
  std::reverse(r.begin(), r.end());
  return r;
}

void addST(const uint32_t n, std::vector<uint32_t>& stree, uint32_t pos, uint32_t value) {
  pos += n;
  stree[pos] += value;
  for (; pos > 1; pos >>= 1) {
    stree[pos>>1] = stree[pos] + stree[pos^1];
  }
}

void addST(const uint32_t n, std::vector<uint32_t>& stree, uint32_t pos) {
  addST(n, stree, pos, 1);
}

uint64_t sumST(const uint32_t n, const std::vector<uint32_t>& stree, uint32_t l, uint32_t r) {
  assert(l <= r);
  uint64_t res = 0;
  for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
    if (l&1) res += stree[l++];
    if (r&1) res += stree[--r];
  }
  return res;
}

bool allocateMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint32_t max_size) {
  if (matrix.size() >= n1 && matrix[0].size() >= n2)
    return false;
  uint32_t new_size1 = std::min(max_size, (uint32_t)(n1 * 1.5));
  uint32_t new_size2 = std::min(max_size, (uint32_t)(n2 * 1.5));
  matrix.clear();
  matrix.resize(new_size1, std::vector<uint64_t>(new_size2));
  return true;
}

bool allocateMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint32_t max_size, uint64_t value) {
  bool result = allocateMatrix(matrix, n1, n2, max_size);
  fillMatrix(matrix, n1, n2, value);
  return result;
}

void fillMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint64_t value) {
  for (uint32_t i = 0; i < n1; i++) {
    for (uint32_t j = 0; j < n2; j++) {
      matrix[i][j] = value;
    }
  }
}
