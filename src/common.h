#pragma once

#include <algorithm>
#include <numeric>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

constexpr double EPS = 1e-8;
constexpr double PI = 3.14159265358979323846;

constexpr uint64_t ONE64 = static_cast<uint64_t>(1);
constexpr uint64_t NOT_SET64 = uint64_t(-1);
constexpr uint32_t NOT_SET32 = uint32_t(-1);
constexpr uint32_t NOT_SET16 = uint16_t(-1);
constexpr uint64_t MAX_64 = std::numeric_limits<uint64_t>::max();
constexpr uint64_t MAX_32 = std::numeric_limits<uint32_t>::max();

template<class T>
inline T Abs(const T& t) {
  if (t > 0)
    return t;
  return -t;
}

template<class T>
inline T AbsDiff(const T& t1, const T& t2) {
  if (t1 >= t2)
    return t1 - t2;
  return t2 - t1;
}

template<class T>
inline T Sgn(const T& t) {
  if (t > 0)
    return 1;
  if (t < 0)
    return -1;
  return 0;
}

template<class T>
inline T Sqr2(const T& t) {
  return ((t) * (t));
}

template <typename T>
std::string to_string(const T& n) {
  std::ostringstream ss;
  ss << n;
  return ss.str();
}

template <typename T>
std::string to_string(const std::vector<T>& vec, const std::string& separator) {
  std::string desc = "";

  for (auto p : vec) {
    if (desc.length() > 0)
      desc += separator;

    desc += to_string(p);
  }

  return desc;
}

template <typename T>
std::string to_string(const std::vector<T>& vec) {
  return to_string(vec, " ");
}

inline int to_int(const std::string& s) {
  int n;
  std::istringstream(s) >> n;
  return n;
}

inline size_t to_sizet(const std::string& s) {
  size_t n;
  std::istringstream(s) >> n;
  return n;
}

inline uint64_t to_uint64(const std::string& s) {
  uint64_t n;
  std::istringstream(s) >> n;
  return n;
}

inline double to_double(const std::string& s) {
  double n;
  std::istringstream(s) >> n;
  return n;
}

inline bool ends_with(const std::string& s, const std::string& postfix) {
  if (postfix.size() > s.size())
    return false;
  return std::equal(postfix.rbegin(), postfix.rend(), s.rbegin());
}

inline bool starts_with(const std::string& s, const std::string& prefix) {
  if (prefix.size() > s.size())
    return false;
  return std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::vector<std::string> SplitNotNull(const std::string& s, const std::string& c);
// std::vector<int> SplitNotNullInt(const std::string& s, const std::string& c);

struct Rand {
  //static size_t setSeed();
  //static size_t setSeed(size_t seed);

  //static double nextDouble();
  static double nextDouble(std::mt19937& rng);
  //static bool check(double probability);
  static bool check(double probability, std::mt19937& rng);
  //static int next();
  //static int next(int bound);
  static int next(int bound, std::mt19937& rng);
  //static int next(int lower, int upper);
  template<typename RandomIt>
  static void shuffle(RandomIt first, RandomIt last, std::mt19937& rng);
  //static std::vector<uint32_t> permutation(uint32_t n);
  static std::vector<uint32_t> permutation(uint32_t n, std::mt19937& rng);
};

template<typename RandomIt>
void Rand::shuffle(RandomIt first, RandomIt last, std::mt19937& rng) {
  typename std::iterator_traits<RandomIt>::difference_type i, n;
  n = last - first;
  for (i = n-1; i > 0; --i) {
    std::swap(first[i], first[Rand::next(i + 1, rng)]);
  }
}

std::vector<uint32_t> identity(uint32_t n);
std::vector<uint32_t> reverse(const std::vector<uint32_t>& p);
std::vector<uint32_t> inverse(const std::vector<uint32_t>& p);

template <typename T>
void sort_unique(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());  
}

template <typename T>
double average(const std::vector<T>& vec) {
  double res = std::accumulate(vec.begin(), vec.end(), 0.0);
  if (!vec.empty())
    res /= double(vec.size());
  return res;
}

template <typename T>
double median(const std::vector<T>& vec) {
  if (vec.empty())
    return 0.0;
  if (vec.size() == 1)
    return vec[0];
  if (vec.size() == 2)
    return (vec[0] + vec[1]) / 2.0;

  std::vector<T> tmp(vec.begin(), vec.end());
  std::sort(tmp.begin(), tmp.end());

  size_t sz = tmp.size();
  if (sz % 2 == 0)
    return (tmp[sz / 2 - 1] + tmp[sz / 2]) / 2.0;
  else
    return tmp[sz / 2];
}

template <typename T>
T Max(const std::vector<T>& v) {
  return *std::max_element(v.begin(), v.end());
}

template <typename T>
T Min(const std::vector<T>& v) {
  return *std::min_element(v.begin(), v.end());
}

double Sum(const std::vector<double>& v);
double Percentile(const std::vector<double>& v, int value);

int Compare(double numberA, double numberB);
bool Equal(double a, double b);
bool GreaterOrEqual(double numberA, double numberB);
bool Greater(double numberA, double numberB);
bool LessOrEqual(double numberA, double numberB);
bool Less(double numberA, double numberB);

/// Segment Tree: increment value at position pos
void addST(const uint32_t n, std::vector<uint32_t>& stree, uint32_t pos);
void addST(const uint32_t n, std::vector<uint32_t>& stree, uint32_t pos, uint32_t value);
/// Segment Tree: sum on interval [l, r)
uint64_t sumST(const uint32_t n, const std::vector<uint32_t>& stree, uint32_t l, uint32_t r);

bool allocateMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint32_t max_size);
bool allocateMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint32_t max_size, uint64_t value);
void fillMatrix(std::vector<std::vector<uint64_t>>& matrix, uint32_t n1, uint32_t n2, uint64_t value);
