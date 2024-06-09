#include "common.h"
#include "logging.h"
#include "cmd_options.h"
#include "bp.h"

#include <csignal>
#include <functional>
#include <iostream>
#include <fstream>
#include <exception>
#include <map>
#include <set>
#include <unordered_set>

volatile sig_atomic_t TLE = 0;
std::function<void(int)> shutdown_handler;
void signal_handler(int signal) { shutdown_handler(signal); }

CMDOptions* Options;

void prepareCMDOptions(CMDOptions& options) {
  options.SetUsageMessage("Usage: pace [options]");

  options.AddAllowedOption("-verbose", "0", "Output debug info");
  options.AddAllowedOption("-seed", "0", "Random seed");

#if defined(HEURISTIC)
  options.AddAllowedOption("-graceful", "1", "Graceful exit");
  options.AddAllowedOption("-confidence", "0", "Minimum confidence to output a solution");
  options.AddAllowedOption("-time-limit", "290", "Time limit in seconds; 0 means no limit");
#elif defined(EXACT)
  options.AddAllowedOption("-graceful", "1", "Graceful exit");
  options.AddAllowedOption("-confidence", "30", "Minimum confidence to output a solution");
  options.AddAllowedOption("-time-limit", "600", "Time limit in seconds; 0 means no limit");
#elif defined(CUTWIDTH)
  options.AddAllowedOption("-graceful", "1", "Graceful exit");
  options.AddAllowedOption("-confidence", "0", "Minimum confidence to output a solution");
  options.AddAllowedOption("-time-limit", "60", "Time limit in seconds; 0 means no limit");
#else
  options.AddAllowedOption("-graceful", "0", "Graceful exit");
  options.AddAllowedOption("-confidence", "0", "Minimum confidence to output a solution");
  options.AddAllowedOption("-time-limit", "0", "Time limit in seconds; 0 means no limit");
#endif

  options.AddAllowedOption("-use-misc", "1", "Whether to apply misc heuristic reordering");
  options.AddAllowedOption("-bf-limit", "10", "The maximum number of vertices to brute-force");

  options.AddAllowedOption("-use-bp", "1", "Whether to apply BP reordering");
  options.AddAllowedOption("-bp-iters", "-1", "The maximum number of BP repetitions (with different seeds)");
  options.AddAllowedOption("-main-bp-iters", "1", "The maximum number of simple-BP repetitions");
  options.AddAllowedOption("-part-bp-iters", "10", "The maximum number of partial-BP repetitions");
  options.AddAllowedOption("-part-bp-size", "768", "Maximum size of partial-BP interval");

  options.AddAllowedOption("-merge-twins", "1", "Merge twin documents");
  options.AddAllowedOption("-split-depth", "18", "The maximum number of recursive splits for BP");
  options.AddAllowedOption("-split-attempts", "1", "The maximum attempts to perform a split");
  options.AddAllowedOption("-iterations-per-split", "40", "The maximum number of BP iterations per split");
  options.AddAllowedOption("-skip-probability", "0.01", "The probability to skip a move");
  options.AddAllowedOption("-backtrack", "1", "Whether to use backtracking for BP reordering");
  options.AddAllowedOption("-use-actual-gains", "false", "Use actual move gains for swapping");

  options.AddAllowedOption("-max-docs-lb", "6144", "The maximum number of documents for computing exact move gains");
  options.AddAllowedOption("-leaf-interval", "18", "The maximum size of leaf intervals to apply full reordering");
  options.AddAllowedOption("-opt-interval-size", "16", "The size of intervals for post-processing of reordered documents");
  options.AddAllowedOption("-opt-interval-iter", "15", "The size of intervals for post-processing of reordered documents");
  options.AddAllowedOption("-post-tune-all", "true", "Post-processing fine-tuning of the complete solution");
  options.AddAllowedOption("-post-tune-int", "true", "Post-processing fine-tuning of the intervals");

  options.AddAllowedOption("-order-alg", "local-median", "Algorithm for ordering at lowest level");
  options.AddAllowedValue("-order-alg", "input");
  options.AddAllowedValue("-order-alg", "median");
  options.AddAllowedValue("-order-alg", "average");
  options.AddAllowedValue("-order-alg", "local-median");
  options.AddAllowedValue("-order-alg", "local-average");
  options.AddAllowedOption("-split-alg", "random", "Algorithm for bisecting at each level");
  options.AddAllowedValue("-split-alg", "input");
  options.AddAllowedValue("-split-alg", "random");
  options.AddAllowedValue("-split-alg", "median");
  options.AddAllowedValue("-split-alg", "average");

  options.AddAllowedOption("-i", "", "Input graph file name");
  options.AddAllowedOption("-o", "", "Output graph file name");
}

struct Result {
  std::string name;
  std::vector<uint32_t> order;
  uint64_t numCrossings = NOT_SET64;
};

struct Graph {
  uint32_t nA; // fixed part
  uint32_t nB; // free part
  uint32_t m;
  uint64_t optCrossings = NOT_SET64; // from the solution file
  std::vector<std::pair<uint32_t, uint32_t>> edges; // u < v, 0-based indices
  Result bestResult;
  std::vector<uint64_t> finalResults;
  uint64_t lowerBound = NOT_SET64;
  uint32_t maxDegreeB;

  std::vector<std::vector<uint32_t>> b2a; // B side to A

  uint32_t origNA;
  uint32_t origNB;
  std::vector<uint32_t> origIndexB;
  std::vector<uint32_t> droppedB;

private:
  mutable std::vector<uint32_t> stree;
  const size_t verbose;

public:
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(size_t verbose): nA(0), nB(0), m(0), verbose(verbose) {}

  void init();
  void simplifyA();
  void simplifyB();
 
  void read(std::istream& in);
  void read(const std::string& filename);

  void printBest(std::ostream& out, const size_t minConfidence) const;
  void printBest(const std::string& filename, const size_t minConfidence) const;

  bool verifyResult(const std::vector<uint32_t>& order) const;

  bool verifyResult(const Result& res) const {
    return verifyResult(res.order);
  }

  uint64_t countCrossings(const std::vector<uint32_t>& order, bool verify=true) const;

  uint64_t countCrossings(const Result& res) const {
    return countCrossings(res.order);
  }

  void setResultIfBetter(const Result& res);

  void setLowerBound(uint64_t lb) {
    if (lb == NOT_SET64)
      return;
    if (lowerBound == NOT_SET64 || lowerBound < lb)
      lowerBound = lb;
  }

  size_t computeConfidence() const;

  void printConfidenceInfo();
};

void Graph::simplifyA() {
  std::vector<uint32_t> degreeA(nA, 0);
  for (const auto& [u, v] : edges) {
    degreeA[u]++;
  }
  size_t numIsolatedA = 0;
  for (size_t i = 0; i < nA; i++) {
    if (degreeA[i] == 0)
      numIsolatedA++;
  }
  if (numIsolatedA == 0) 
    return;

  // Get rid of isolated vertices on side A: drop the vertices and re-number edges
  CHECK(numIsolatedA < nA);

  std::vector<uint32_t> newIndexA;
  uint32_t curIdx = 0;
  for (size_t i = 0; i < nA; i++) {
    if (degreeA[i] == 0)
      newIndexA.push_back(NOT_SET32);
    else
      newIndexA.push_back(curIdx++);
  }
  for (size_t i = 0; i < edges.size(); i++) {
    uint32_t u = edges[i].first;
    CHECK(degreeA[u] > 0 && newIndexA[u] != NOT_SET32);
    edges[i].first = newIndexA[u];
  }

  nA -= numIsolatedA;
  CHECK(curIdx == nA);
  LOG_IF(verbose, "removed %d isolated A vertices; new nA = %d", numIsolatedA, nA);
}

void Graph::simplifyB() {
  std::vector<uint32_t> degreeB(nB, 0);
  for (const auto& [u, v] : edges) {
    degreeB[v]++;
  }
  size_t numIsolatedB = 0;
  for (size_t i = 0; i < nB; i++) {
    if (degreeB[i] == 0)
      numIsolatedB++;
  }
  if (numIsolatedB == 0) 
    return;

  // Get rid of isolated vertices on side B: keep the dropped indices and re-number edges
  CHECK(numIsolatedB < nB);
  origNB = nB;
  origIndexB.resize(nB, NOT_SET32);

  std::vector<uint32_t> newIndexB;
  uint32_t curIdx = 0;
  for (size_t i = 0; i < nB; i++) {
    if (degreeB[i] == 0) {
      newIndexB.push_back(NOT_SET32);
      droppedB.push_back(i);
    } else {
      origIndexB[curIdx] = i;
      newIndexB.push_back(curIdx++);
    }
  }
  for (size_t i = 0; i < edges.size(); i++) {
    uint32_t u = edges[i].second;
    CHECK(degreeB[u] > 0 && newIndexB[u] != NOT_SET32);
    edges[i].second = newIndexB[u];
  }

  nB -= numIsolatedB;
  CHECK(curIdx == nB);
  LOG_IF(verbose, "removed %d isolated B vertices; new nB = %d", numIsolatedB, nB);
}

void Graph::init() {
  CHECK(nA > 0 && nB > 0 && m > 0);

  simplifyA();
  simplifyB();

  // check the degrees
  std::vector<uint32_t> degreeA(nA, 0);
  std::vector<uint32_t> degreeB(nB, 0);
  std::vector<uint32_t> singleAdjA(nA, NOT_SET32);
  for (const auto& [u, v] : edges) {
    degreeA[u]++;
    degreeB[v]++;
    singleAdjA[u] = v;
  }
  LOG_IF(verbose, "%'lld <= degree_A <= %'lld;  %'lld <= degree_B <= %'lld;", Min(degreeA), Max(degreeA), Min(degreeB), Max(degreeB));
  maxDegreeB = Max(degreeB);

  size_t numIsolatedA = 0;
  size_t numDegreeOneA = 0;
  for (size_t i = 0; i < nA; i++) {
    if (degreeA[i] == 0)
      numIsolatedA++;
    if (degreeA[i] == 1)
      numDegreeOneA++;
  }
  (void)numIsolatedA;
  size_t numIsolatedB = 0;
  size_t numDegreeOneB = 0;
  for (size_t i = 0; i < nB; i++) {
    if (degreeB[i] == 0)
      numIsolatedB++;
    if (degreeB[i] == 1)
      numDegreeOneB++;
  }
  (void)numIsolatedB;
  CHECK(numIsolatedA == 0 && numIsolatedB == 0, "smth wrong with simplification");
  LOG_IF(verbose, "degree-1 B vertices: %d; degree-1 A vertices: %d", numDegreeOneB, numDegreeOneA);

  // build the data structures
  stree.resize(2 * nA);

  b2a.resize(nB);
  for (size_t i = 0; i < nB; i++) {
    b2a[i].reserve(degreeB[i]);
  }

  for (const auto& [u, v] : edges) {
    CHECK(0 <= u && u < nA);
    CHECK(0 <= v && v < nB);
    b2a[v].push_back(u);
  }

  size_t numMultiEdges = 0;
  for (size_t i = 0; i < nB; i++) {
    std::sort(b2a[i].begin(), b2a[i].end());
    uint32_t prevU = NOT_SET32;
    for (uint32_t u : b2a[i]) {
      CHECK(prevU == NOT_SET32 || prevU <= u);
      if (prevU == u)
        numMultiEdges++;
      prevU = u;
    }
  }
  LOG_IF(verbose && numMultiEdges, "multi-edges: %d", numMultiEdges);
  // CHECK(numMultiEdges == 0);
}

void Graph::setResultIfBetter(const Result& res) {
  CHECK(res.numCrossings != NOT_SET64);
  std::string showName = "'" + res.name + "'";
  if (bestResult.numCrossings == NOT_SET64 || bestResult.numCrossings >= res.numCrossings) {
    bestResult = res;
    LOG_IF(verbose, "solution %-18s has \033[95m%'lld\033[0m crossings", showName.c_str(), res.numCrossings);
  } else {
    LOG_IF(verbose, "solution %-18s has %'lld crossings", showName.c_str(), res.numCrossings);
  }

  if (ends_with(res.name, "final")) {
    finalResults.push_back(res.numCrossings);
  }
}

size_t Graph::computeConfidence() const {
  if (bestResult.numCrossings <= lowerBound) {
    return 100;
  }

  size_t numOptimal = 0;
  for (uint64_t res : finalResults) {
    if (res == bestResult.numCrossings)
      numOptimal++;
  }
  size_t numAll = std::max(finalResults.size(), size_t(10));
  size_t conf = size_t(100.0 * double(numOptimal) / numAll);
  LOG_IF(verbose, "  confidence: %d%% (%d out of %d)", conf, numOptimal, finalResults.size());
  return conf;
}

void Graph::printConfidenceInfo() {
  if (!verbose) {
    return;
  }

  if (bestResult.numCrossings <= lowerBound) {
    LOG_IF(verbose, "  confidence: %.1lf%% (reached lower bound)", 100.0);
    return;
  }
  if (finalResults.empty()) {
    LOG_IF(verbose, "  confidence: none");
    return;
  }

  computeConfidence();
}

uint64_t Graph::countCrossings(const std::vector<uint32_t>& order, bool verify) const {
  CHECK(!verify || verifyResult(order));

  const uint32_t n = nA;
  std::fill(stree.begin(), stree.end(), 0);

  uint64_t numCrossFast = 0;
  for (uint32_t v : order) {  
    // check every edge for v
    for (uint32_t u : b2a[v]) {
      numCrossFast += sumST(n, stree, u + 1, n);
    }
    // add the edges
    for (uint32_t u : b2a[v]) {
      addST(n, stree, u);
    }
  }

  if (verify && verbose >= 2 && order.size() <= 1024) {
    // debug-only slow version
    auto idx = inverse(order);
    std::vector<std::pair<uint32_t, uint32_t>> redges;
    for (const auto& [u, v] : edges) {
      redges.push_back({u, idx[v]});
    }
    uint64_t numCrossSlow = 0;
    for (size_t i = 0; i < redges.size(); i++) {
      for (size_t j = i + 1; j < redges.size(); j++) {
        if ((redges[i].first < redges[j].first && redges[i].second > redges[j].second) ||
            (redges[i].first > redges[j].first && redges[i].second < redges[j].second))
          numCrossSlow++;
      }
    }

    if (numCrossSlow != numCrossFast)
      LOG_IF(verbose, "order = %s", to_string(order).c_str());
    CHECK(numCrossSlow == numCrossFast, "numCrossSlow = %d, numCrossFast = %d", numCrossSlow, numCrossFast);
  }

  return numCrossFast;
}

bool Graph::verifyResult(const std::vector<uint32_t>& order) const {
  if (nB != order.size())
    return false;
  std::vector<bool> taken(nB, false);
  for (size_t i = 0; i < nB; i++) {
    if (!(0 <= order[i] && order[i] < nB))
      return false;
    if (taken[order[i]])
      return false;
    taken[order[i]] = true;
  }
  return true;
}

void Graph::read(std::istream& in) {
  std::string line;  
  while (std::getline(in, line)) {
    if (starts_with(line, "c opt")) {
      CHECK(optCrossings == NOT_SET64);
      auto tmp = SplitNotNull(line, " ");
      optCrossings = to_uint64(tmp[2]);
      CHECK(optCrossings >= 0);
    }
    if (line[0] == 'c')
      continue;
    if (line.length() == 0)
      continue;
    if (starts_with(line, "p ocr")) {
      CHECK(nA == 0 && nB == 0 && m == 0);
      auto tmp = SplitNotNull(line, " ");
      nA = to_int(tmp[2]);
      nB = to_int(tmp[3]);
      m = to_int(tmp[4]);
      CHECK(nA > 0 && nB > 0 && m > 0);

      // check if that's a cutwidth
      if (tmp.size() == 6) {
        for (size_t i = 0; i < nA + nB; i++) {
          std::getline(in, line);
        }
      }
      continue;
    }
    CHECK(nA > 0 && nB > 0 && m > 0);
    auto tmp = SplitNotNull(line, " ");
    uint32_t u = to_int(tmp[0]);
    uint32_t v = to_int(tmp[1]);
    CHECK(1 <= u && u <= nA);
    CHECK(nA + 1 <= v && v <= nA + nB);
    edges.push_back({u - 1, v - 1 - nA});
  }
  CHECK(edges.size() == m);

  LOG_IF(verbose, "parsed graph with nA = %'lld, nB = %'lld, and |E| = %'lld", nA, nB, edges.size());
  origNA = nA;
  origNB = nB;
}

void Graph::read(const std::string& filename) {
  if (filename.empty()) {
    read(std::cin);
    return;
  }

  std::ifstream fileStream;
  fileStream.open(filename.c_str(), std::ios::in);
  if (!fileStream)
    ERROR("input file '" + filename + "' doesn't exist");
  LOG_IF(verbose, "processing '%s'", filename.c_str());
  read(fileStream);
  fileStream.close();
}

void Graph::printBest(std::ostream& out, const size_t minConfidence) const {
  if (minConfidence > 0) {
    const size_t confidence = computeConfidence();
    if (confidence < minConfidence) {
      LOG_IF(verbose, "skipping printing a solution due to low confidence (%d < %d)", confidence, minConfidence);
      exit(1);
    }
  }

  if (!verifyResult(bestResult)) {
    LOG_IF(verbose, "printing backup solution");
    for (size_t i = 0; i < origNB; i++) {
      out << i + origNA + 1 << "\n";
    }
    return;
  }

  LOG_IF(verbose, "printing solution '%s'", bestResult.name.c_str());
  if (droppedB.empty()) {
    for (size_t i = 0; i < bestResult.order.size(); i++) {
      out << bestResult.order[i] + origNA + 1 << "\n";
    }
  } else {
    for (size_t i = 0; i < bestResult.order.size(); i++) {
      uint32_t v = bestResult.order[i];
      out << origIndexB[v] + origNA + 1 << "\n";
    }
    for (size_t i = 0; i < droppedB.size(); i++) {
      out << droppedB[i] + origNA + 1 << "\n";
    }
  }
}

void Graph::printBest(const std::string& filename, const size_t minConfidence) const {
  if (filename.empty()) {
    printBest(std::cout, minConfidence);
    return;
  }

  std::ofstream fileStream;
  fileStream.open(filename.c_str(), std::ofstream::out);
  printBest(fileStream, minConfidence);
  fileStream.close();
}

void applyIdentity(Graph& graph, size_t verbose) {
  Result res;
  res.name = "identity";
  res.order = identity(graph.nB);
  res.numCrossings = graph.countCrossings(res);
  graph.setResultIfBetter(res);
}

void applyReverse(Graph& graph, size_t verbose) {
  Result res;
  res.name = "reverse";
  res.order = reverse(identity(graph.nB));
  res.numCrossings = graph.countCrossings(res);
  graph.setResultIfBetter(res);
}

void applyBruteForce(Graph& graph, size_t verbose) {
  LOG_IF(verbose, "running brute-force for %d vertices ...", graph.nB);
  std::vector<uint32_t> p = identity(graph.nB);
  size_t iter = 0;
  do {
    Result res;
    res.name = "bf_" + std::to_string(iter);
    res.order = p;
    res.numCrossings = graph.countCrossings(res);
    graph.setResultIfBetter(res);
    iter++;
  } while (std::next_permutation(p.begin(), p.end()));
}

void applyMedian(Graph& graph, size_t verbose) {
  std::vector<double> med(graph.nB);
  std::vector<double> avg(graph.nB);
  for (size_t i = 0; i < graph.nB; i++) {
    med[i] = median(graph.b2a[i]);
    avg[i] = average(graph.b2a[i]);
  }
  std::vector<uint32_t> order = identity(graph.nB);
  std::stable_sort(order.begin(), order.end(), [&](uint32_t l, uint32_t r) {
    return std::make_tuple(med[l], avg[l], l) < std::make_tuple(med[r], avg[r], r);
  });

  Result res;
  res.name = "median";
  res.order = order;
  res.numCrossings = graph.countCrossings(res);
  graph.setResultIfBetter(res);
  CHECK(graph.optCrossings == NOT_SET64 || res.numCrossings <= 3 * graph.optCrossings);
}

void applyAverage(Graph& graph, size_t verbose) {
  std::vector<double> med(graph.nB);
  std::vector<double> avg(graph.nB);
  for (size_t i = 0; i < graph.nB; i++) {
    med[i] = median(graph.b2a[i]);
    avg[i] = average(graph.b2a[i]);
  }
  std::vector<uint32_t> order = identity(graph.nB);
  std::stable_sort(order.begin(), order.end(), [&](uint32_t l, uint32_t r) {
    return std::make_tuple(avg[l], med[l], l) < std::make_tuple(avg[r], med[r], r);
  });

  Result res;
  res.name = "average";
  res.order = order;
  res.numCrossings = graph.countCrossings(res);
  graph.setResultIfBetter(res);
}

void printResult(const std::vector<Document *>& documentsPtr, uint32_t numDocs, const std::string& filename) {
  std::ofstream out;
  out.open(filename.c_str(), std::ios::out);
  size_t numEdges = 0;
  std::set<uint32_t> utils;
  for (size_t i = 0; i < numDocs; i++) {
    Document* doc = documentsPtr[i];
    for (const auto [u, _] : doc->getEdges()) {
      utils.insert(u);
      numEdges++;
    }
  }
  size_t numUtilities = 0;
  std::unordered_map<uint32_t, uint32_t> u2i;
  for (uint32_t u : utils) {
    u2i[u] = numUtilities++;
  }
  out << "p ocr " << numUtilities << " " << numDocs << " " << numEdges << "\n";
  for (size_t i = 0; i < numDocs; i++) {
    Document* doc = documentsPtr[i];
    for (const auto [u, _] : doc->getEdges()) {
      out << u2i[u] + 1 << " " << numUtilities + i + 1 << "\n";
    }
  }
  out.close();
}

void findDuplicates(CMDOptions& options, Graph& graph, std::vector<std::vector<uint32_t>>& doc2ids) {
  const size_t verbose = options.getInt("-verbose");

  doc2ids.reserve(graph.nB);
  std::map<std::vector<uint32_t>, uint32_t> edge2head;
  size_t numDuplicates = 0;
  for (uint32_t i = 0; i < graph.nB; i++) {
    const std::vector<uint32_t>& edges = graph.b2a[i];
    CHECK(!edges.empty());

    auto it = edge2head.find(edges);
    if (it == edge2head.end()) {
      // Found a new doc
      uint32_t head = doc2ids.size();
      doc2ids.push_back({i});
      edge2head[edges] = head;
    } else {
      // Found a duplicate
      uint32_t head = it->second;
      doc2ids[head].push_back(i);
      numDuplicates++;
    }
  }

  CHECK(numDuplicates + doc2ids.size() == graph.nB);
  LOG_IF(verbose, "  merged %3d duplicates; remaining docs: %d", numDuplicates, doc2ids.size());
}

void findTwins(CMDOptions& options, Graph& graph, std::vector<std::vector<uint32_t>>& doc2ids) {
  const size_t verbose = options.getInt("-verbose");

  std::vector<std::vector<uint32_t>> duplicates = doc2ids;
  doc2ids.clear();
  std::vector<uint32_t> heads;

  // Compute degrees
  std::vector<uint32_t> degreeU(graph.nA, 0);
  for (uint32_t i = 0; i < duplicates.size(); i++) {
    const std::vector<uint32_t>& ids = duplicates[i];
    CHECK(!ids.empty());

    const std::vector<uint32_t>& edges = graph.b2a[ids[0]];
    for (uint32_t u : edges) {
      CHECK(u < degreeU.size());
      degreeU[u]++;
    }
  }

  // Do the edges contain the utility?
  auto hasEdge = [](const std::vector<uint32_t>& edges, uint32_t utility) {
    for (uint32_t u : edges) {
      if (u == utility)
        return true;
    }
    return false;
  };

  // Are the two vertyices twins?
  auto isTwin = [&](uint32_t i, uint32_t j) {
    CHECK(i < graph.nB);
    CHECK(j < graph.nB);
    const std::vector<uint32_t>& edges1 = graph.b2a[i];
    const std::vector<uint32_t>& edges2 = graph.b2a[j];
    if (edges1.size() != edges2.size())
      return false;

    for (uint32_t k = 0; k < edges1.size(); k++) {
      const uint32_t u1 = edges1[k];
      const uint32_t u2 = edges2[k];

      uint32_t degU1 = degreeU[u1] - 1;
      // reduce degree from doc1 and doc2
      CHECK(hasEdge(edges1, u1) && degreeU[u1] >= 1);
      if (hasEdge(edges2, u1)) {
        CHECK(degU1 >= 1);
        degU1--;
      }

      uint32_t degU2 = degreeU[u2] - 1;
      CHECK(hasEdge(edges2, u2) && degreeU[u2] >= 1);
      if (hasEdge(edges1, u2)) {
        CHECK(degU2 >= 1);
        degU2--;
      }

      if (u1 == u2 || (AbsDiff(u1, u2) == 1 && degU1 == 0 && degU2 == 0)) {
        // pass
      } else {
        return false;
      }
    }
    return true;
  };

  for (uint32_t i = 0; i < duplicates.size(); i++) {
    const std::vector<uint32_t>& ids = duplicates[i];
    CHECK(!ids.empty());

    uint32_t twinIdx = NOT_SET32;
    // split heads by degree;
    if (heads.size() <= 4096) {
      for (uint32_t j = 0; j < heads.size(); j++) {
        // maybe check weighted twins
        if (isTwin(heads[j], ids[0])) {
          twinIdx = j;
          break;
        }
      }
    }

    if (twinIdx == NOT_SET32) {
      doc2ids.push_back(ids);
      heads.push_back(ids[0]);
    } else {
      doc2ids[twinIdx].insert(doc2ids[twinIdx].end(), ids.begin(), ids.end());
    }
  }

  CHECK(duplicates.size() >= doc2ids.size());
  size_t numTwins = duplicates.size() - doc2ids.size();
  LOG_IF(verbose, "  merged %3d      twins; remaining docs: %d", numTwins, doc2ids.size());
}

uint64_t sortTwins(CMDOptions& options, Graph& graph, std::vector<uint32_t>& ids) {
  const size_t verbose = options.getInt("-verbose");
  CHECK(ids.size() > 1);
  std::unordered_map<uint32_t, double> med;
  med.reserve(ids.size());
  std::unordered_map<uint32_t, double> avg;
  avg.reserve(ids.size());
  for (size_t i = 0; i < ids.size(); i++) {
    const uint32_t id = ids[i];
    med[id] = median(graph.b2a[id]);
    avg[id] = average(graph.b2a[id]);
  }
  // median
  uint64_t medCrossings = 0;
  if (verbose) {
    std::sort(ids.begin(), ids.end(), [&](uint32_t l, uint32_t r) {
      return std::make_tuple(med[l], avg[l], l) < std::make_tuple(med[r], avg[r], r);
    });
    medCrossings = graph.countCrossings(ids, false);
    (void)medCrossings;
  }
  // average
  std::sort(ids.begin(), ids.end(), [&](uint32_t l, uint32_t r) {
    return std::make_tuple(avg[l], med[l], l) < std::make_tuple(avg[r], med[r], r);
  });
  uint64_t avgCrossings = graph.countCrossings(ids, false);
  if (verbose) {
    CHECK(avgCrossings == medCrossings);
  }

  return avgCrossings;
}

void initDocuments(CMDOptions& options, Graph& graph, std::vector<Document>& documents) {
  const size_t verbose = options.getInt("-verbose");

  LOG_IF(verbose, "kernelization:");

  // Find unique documents
  std::vector<std::vector<uint32_t>> doc2ids;
  findDuplicates(options, graph, doc2ids);
  if (options.getBool("-merge-twins")) {
    findTwins(options, graph, doc2ids);
  }

  const uint32_t numDocs = doc2ids.size();
  documents.resize(numDocs);
  std::map<uint32_t, uint32_t> edgeCounts;
  uint32_t maxUtility = 0;
  uint32_t numIds = 0;

  // Init documents and edges
  for (uint32_t i = 0; i < numDocs; i++) {
    const auto& ids = doc2ids[i];
    CHECK(!ids.empty());
    numIds += ids.size();
    Document& doc = documents[i];
    doc.Index = i;
    if (ids.size() == 1) {
      doc.init(ids);
    } else {
      std::vector<uint32_t> sortedIds = ids;
      uint64_t crossings = sortTwins(options, graph, sortedIds);
      doc.setSelfCrossings(crossings);
      doc.init(sortedIds);
    }

    // Collect doc edges
    edgeCounts.clear();
    for (uint32_t j = 0; j < ids.size(); j++) {
      for (uint32_t u : graph.b2a[ids[0]]) {
        edgeCounts[u] += 1;
        maxUtility = std::max(maxUtility, u);
      }
    }
    std::vector<EdgeTy> edges;
    edges.reserve(edgeCounts.size());
    for (const auto [u, w] : edgeCounts) {
      edges.emplace_back(u, w);
    }
    doc.setOrigEdges(edges);
  }
  (void)numIds;
  CHECK(graph.nB == numIds, "graph.nB = %d; numIds = %d", graph.nB, numIds);

  ///////////////////////// Extra Kernelization ////////////////////////////////
  // Merge edges to consecutive degree-1 A vertices
  std::vector<uint32_t> utililtyDegree(maxUtility + 1, 0);
  for (uint32_t i = 0; i < documents.size(); i++) {
    const Document& doc = documents[i];
    for (const auto [u, _] : doc.getEdges()) {
      utililtyDegree[u]++;
    }
  }
  size_t numMergedDegree1A = 0;
  for (uint32_t i = 0; i < documents.size(); i++) {
    Document& doc = documents[i];
    std::vector<EdgeTy> newEdges;
    newEdges.reserve(doc.getEdges().size());
    uint32_t prevU = NOT_SET32;
    for (const auto [u, w] : doc.getEdges()) {
      if (!newEdges.empty() && utililtyDegree[newEdges.back().utility] == 1 && utililtyDegree[u] == 1 && (prevU != NOT_SET32 && prevU + 1 == u)) {
        newEdges.back().weight += w;
        numMergedDegree1A++;
        prevU = u;
        continue;
      }
      newEdges.emplace_back(u, w);
      prevU = u;
    }
    if (newEdges.size() != doc.getEdges().size()) {
      newEdges.shrink_to_fit();
      doc.setOrigEdges(newEdges);
    }
  }
  LOG_IF(verbose && numMergedDegree1A > 0, "  merged %d edges to consecutive degree-1 A vertices", numMergedDegree1A);
  //////////////////////////////////////////////////////////////////////////////

  // Print some stats
  uint32_t numEdges = 0;
  uint32_t numEdgesW = 0;
  std::unordered_set<uint32_t> uniqueUtils;
  for (uint32_t idx = 0; idx < documents.size(); idx++) {
    Document &doc = documents[idx];
    CHECK(doc.Index == idx);
    for (const auto [u, w] : doc.getEdges()) {
      uniqueUtils.insert(u);
      numEdges += 1;
      numEdgesW += w;
    }
  }

  LOG_IF(verbose >= 1, "initialized %'lld documents and %'lld utilities with %'lld edges (%'lld weighted)", numDocs, uniqueUtils.size(), numEdges, numEdgesW);
}

template <class T>
uint64_t mainBP(
    Graph& graph,
    CMDOptions& options,
    BalancedPartitioningConfig& config,
    BalancedPartitioning<T>& alg,
    std::vector<Document *>& documents,
    std::function<void(const std::string& name)> progressFunc,
    std::function<bool()> convergedFunc,
    const OptLevelT OptLevel) {
  const size_t verbose = options.getInt("-verbose");
  const size_t mainBPIters = options.getInt("-main-bp-iters");
  CHECK(mainBPIters >= 1);

  // Run the reordering
  alg.run(documents);

  if (mainBPIters > 1 && OptLevel >= OptLevelT::O3) {
    std::vector<std::pair<uint32_t, uint32_t>> bestBuckets; // (id,bucket)
    collectBuckets(documents, bestBuckets);
    auto order = extractOrderedIds(documents);
    uint64_t bestCrossings = graph.countCrossings(order, false);
    LOG_IF(verbose, "Orig BestCrossings: %'lld", bestCrossings);

    std::mt19937 RNG;
    const size_t seed0 = alg.getSeed();
    for (size_t iter = 1; iter < mainBPIters; iter++) {
      const size_t iterSeed = seed0 * 24524831 + iter;
      RNG.seed(iterSeed);
      Rand::shuffle(documents.begin(), documents.end(), RNG);
      alg.setSeed(iterSeed);
      alg.run(documents);
      alg.setSeed(seed0);

      auto order = extractOrderedIds(documents);
      uint64_t curCrossings = graph.countCrossings(order, false);
      LOG_IF(verbose, "Crossings at iteration %d: %'lld", iter, curCrossings);
      if (curCrossings < bestCrossings) {
        bestCrossings = curCrossings;
        collectBuckets(documents, bestBuckets);
      }
    }
    assignBuckets(documents, bestBuckets);
  }

  CHECK(bucketsCorrect(documents.begin(), documents.end()));

  // Save the result
  progressFunc("bp_" + std::to_string(alg.getSeed()));

  // Apply post-processing improving the solution
  if (!convergedFunc() && config.PostTuneAll) {
    if (config.PostTuneInt) {
      alg.tuneSolution(documents, progressFunc, convergedFunc, OptLevel, alg.getVerbose());
    } else if (OptLevel >= OptLevelT::O3) {
      alg.tuneLargeSolution(documents, progressFunc, convergedFunc, OptLevelT::O1, alg.getVerbose());
    }
  }

  auto order = extractOrderedIds(documents);
  return graph.countCrossings(order, false);
}

template <class T>
uint64_t partialBP(
    Graph& graph,
    CMDOptions& options,
    const size_t seed,
    BalancedPartitioningConfig& config,
    BalancedPartitioning<T>& alg,
    std::vector<Document *>& documents,
    std::function<void(const std::string& name)> progressFunc,
    std::function<bool()> convergedFunc,
    uint32_t First, uint32_t Last) {
  const size_t verbose = options.getInt("-verbose");
  const uint32_t MaxPartBPSize = options.getInt("-part-bp-size");
  while (Last - First > MaxPartBPSize) {
    Last -= 1;
    First += 1;
  }
  CHECK(First < Last);

  const size_t partBPIters = options.getInt("-part-bp-iters");
  CHECK(partBPIters > 0);
  auto beginTime = std::chrono::steady_clock::now();

  const uint32_t numDocs = documents.size();
  LOG_IF(verbose, "running partialBP [%.2lf .. %.2lf] with %4d docs using %d iterations", double(First) / numDocs, double(Last) / numDocs, Last - First, partBPIters);

  CHECK(bucketsCorrect(documents.begin(), documents.end()));
  const uint32_t MinBucket = documents[First]->Bucket;

  auto origOrder = extractOrderedIds(documents);
  uint64_t AllCrossings = graph.countCrossings(origOrder, false);

  std::vector<Document*> subDocs(documents.begin() + First, documents.begin() + Last);
  std::vector<uint32_t> subOrder = extractOrderedIds(subDocs);
  const uint64_t InitSubCrossings = graph.countCrossings(subOrder, false);
  uint64_t BestSubCrossings = InitSubCrossings;
  std::vector<uint64_t> IterObjective = {InitSubCrossings};

  std::mt19937 RNG;
  for (size_t i = 0; i < partBPIters; i++) {
    const size_t iterSeed = seed * 12561 + i;
    RNG.seed(iterSeed);
    Rand::shuffle(subDocs.begin(), subDocs.end(), RNG);

    options.setInt("-verbose", 0);
    alg.setVerbose(0);
    alg.setSeed(iterSeed);
    auto intProgressFunc = [](const std::string& name) {};
    uint64_t NewSubCrossings = mainBP(graph, options, config, alg, subDocs, intProgressFunc, convergedFunc, OptLevelT::O1);
    options.setInt("-verbose", verbose);
    alg.setVerbose(verbose);
    alg.setSeed(seed);

    if (NewSubCrossings < InitSubCrossings)
      LOG_IF(verbose, "  new crossings at iteration %d: %'lld (-%d)", i, AllCrossings - InitSubCrossings + NewSubCrossings, InitSubCrossings - NewSubCrossings);
    else
      LOG_IF(verbose, "  new crossings at iteration %d: %'lld (+%d)", i, AllCrossings - InitSubCrossings + NewSubCrossings, NewSubCrossings - InitSubCrossings);

    if (NewSubCrossings < BestSubCrossings) {
      BestSubCrossings = NewSubCrossings;
      // Set new buckets
      std::sort(documents.begin() + First, documents.begin() + Last, Document::BucketOrder());
      for (uint32_t I = First; I < Last; I++) {
        CHECK(0 <= documents[I]->Bucket && documents[I]->Bucket < subDocs.size());
        documents[I]->Bucket += MinBucket;
      }
      // Save the progress
      progressFunc("bp_" + std::to_string(seed) + "_partialBP-" + std::to_string(i));
    } else {
      // Revert buckets to original
      for (uint32_t I = First; I < Last; I++) {
        documents[I]->Bucket = MinBucket + (I - First);
      }
    }

    auto nowTime = std::chrono::steady_clock::now();
    auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(nowTime - beginTime).count();
    if (durationSec >= 30)
      break;

    IterObjective.push_back(NewSubCrossings);
    if (attemptsConverged(
      IterObjective,
      /* MaxAttempts */ partBPIters,
      /* MaxItersWithoutProgress */ 5
    ))
      break;
  }

  CHECK(bucketsCorrect(documents.begin(), documents.end()));
  CHECK(InitSubCrossings >= BestSubCrossings);
  const uint64_t Delta = InitSubCrossings - BestSubCrossings;
  auto endTime = std::chrono::steady_clock::now();
  auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - beginTime).count();
  if (Delta > 0) {
    LOG_IF(verbose, "  completed partialBP; reduced crossings by %d (%'lld -> %'lld); running time is %d seconds", Delta, AllCrossings, AllCrossings - Delta, durationSec);
    progressFunc("bp_" + std::to_string(seed) + "_partialBP");
  } else {
    LOG_IF(verbose, "  completed partialBP; didn't reduce crossings; running time is %d seconds", durationSec);
  }
  AllCrossings -= Delta;
  return Delta;
}

template <class T>
void createAndApplyBP(
    Graph& graph,
    CMDOptions& options,
    const size_t seed,
    BalancedPartitioningConfig& config,
    BalancedPartitioning<T>& alg,
    std::vector<Document *>& documents) {
  const size_t verbose = options.getInt("-verbose");
  const uint32_t numDocs = documents.size();

  // Init the reordering algorithm
  alg.setVerbose(verbose);
  alg.setSeed(seed);

  // Init the LB data
  alg.initializeLB(documents);

  // Count crossings and remember the solution.
  auto progressFunc = [&](const std::string& name) {
    Result res;
    res.name = name;
    res.order = extractOrderedIds(documents);
    res.numCrossings = graph.countCrossings(res);
    graph.setResultIfBetter(res);
    graph.setLowerBound(alg.getLowerBound());
  };
  // Check alg convergence
  auto convergedFunc = [&graph]() -> bool { return graph.lowerBound >= graph.bestResult.numCrossings; };

  // Run the algorithm
  mainBP(graph, options, config, alg, documents, progressFunc, convergedFunc, OptLevelT::O3);

  const size_t partBPIters = options.getInt("-part-bp-iters");

  uint64_t SavedSubCrossings = 0;
  if (partBPIters > 0 && config.PostTuneInt) {
    const std::vector<std::pair<uint32_t, uint32_t>> subIntervals = {
      {numDocs / 4, 3 * numDocs / 4},
      {0, numDocs / 3},
      {2 * numDocs / 3, numDocs},
      {numDocs / 3, 2 * numDocs / 3}
    };
    for (const auto& [first, last] : subIntervals) {
      SavedSubCrossings += partialBP(graph, options, seed, config, alg, documents, progressFunc, convergedFunc, first, last);
    }

    LOG_IF(verbose, "saved a total of %'lld crossings using partial BP", SavedSubCrossings);
  }

  if (SavedSubCrossings > 0 && !convergedFunc() && config.PostTuneAll) {
    alg.tuneSolution(documents, progressFunc, convergedFunc, OptLevelT::O3, verbose);
  }

  progressFunc("bp_" + std::to_string(seed) + "_final");
}

void applyBP(Graph& graph, CMDOptions& options, size_t seed, BalancedPartitioningConfig& config) {
  // Create and initialize a bipartite graph in which one part is a given
  // set of documents and another part contains the corresponding utility nodes
  std::vector<Document> documents;
  initDocuments(options, graph, documents);
  const uint32_t numDocs = documents.size();

  // Copying the list of documents as it will be modified during the execution
  std::vector<Document *> documentsPtr;
  documentsPtr.reserve(documents.size());
  for (uint32_t idx = 0; idx < documents.size(); idx++) {
    documentsPtr.push_back(&documents[idx]);
  }
  // Drop empty utilities after kernelization
  updateUtilities(documentsPtr.begin(), documentsPtr.end(), true);

  uint32_t numEdges = 0;
  uint32_t maxDocDegree = 0;
  for (uint32_t idx = 0; idx < documents.size(); idx++) {
    numEdges += documents[idx].getEdges().size();
    maxDocDegree = std::max(maxDocDegree, uint32_t(documents[idx].weightedDegree()));
  }

  // Adjust config to avoid excessive runtimes
  adjustBPConfig(options.getInt("-verbose"), config, numDocs, numEdges);

  // Create and initialize LBData
  const bool canUseDense = numDocs <= std::max(config.MaxDocsLB, uint32_t(32));
  const bool canUse8 = uint64_t(maxDocDegree) * uint64_t(maxDocDegree) < uint64_t(std::numeric_limits<uint8_t>::max());

  if (canUse8) {
    auto alg = std::make_unique<BalancedPartitioning<uint8_t>>(config, canUseDense);
    createAndApplyBP<uint8_t>(graph, options, seed, config, *alg, documentsPtr);
    return;
  }

  const bool canUse16 = uint64_t(maxDocDegree) * uint64_t(maxDocDegree) < uint64_t(std::numeric_limits<uint16_t>::max());
  if (canUse16) {
    auto alg = std::make_unique<BalancedPartitioning<uint16_t>>(config, canUseDense);
    createAndApplyBP<uint16_t>(graph, options, seed, config, *alg, documentsPtr);
    return;
  }

  const bool canUse32 = uint64_t(maxDocDegree) * uint64_t(maxDocDegree) < uint64_t(std::numeric_limits<uint32_t>::max());
  if (canUse32) {
    auto alg = std::make_unique<BalancedPartitioning<uint32_t>>(config, canUseDense);
    createAndApplyBP<uint32_t>(graph, options, seed, config, *alg, documentsPtr);
    return;
  }

  auto alg = std::make_unique<BalancedPartitioning<uint64_t>>(config, canUseDense);
  createAndApplyBP<uint64_t>(graph, options, seed, config, *alg, documentsPtr);
}

void adjustBPConfig(int verbose, BalancedPartitioningConfig& config, uint32_t numDocs, uint32_t numEdges) {
  // Size of intervals in post-processing
  const uint64_t MAX_SORT_OPS = ONE64 << 32;
  auto estimatedSortOPS = [&](uint32_t intervalSize) -> uint64_t {
    return (static_cast<uint64_t>(1) << intervalSize) * uint64_t(intervalSize) * uint64_t(numDocs);
  };
  uint32_t intervalSize = config.PostIntervalSize;
  if (intervalSize > numDocs) {
    intervalSize = numDocs;
    config.PostIntervalSize = numDocs;
  }
  if (estimatedSortOPS(intervalSize) > MAX_SORT_OPS) {
    while (intervalSize > 1 && estimatedSortOPS(intervalSize) > MAX_SORT_OPS) {
      intervalSize -= 1;
    }
    LOG_IF(verbose, "reduced interval-size from %d to %d; ops = %'lld", config.PostIntervalSize, intervalSize, estimatedSortOPS(intervalSize));
    config.PostIntervalSize = intervalSize;
  } else {
    LOG_IF(verbose, "using interval-size = %d; ops = %'lld", intervalSize, estimatedSortOPS(intervalSize));
  }

  // Number of documents with exact move gains
  const uint64_t MAX_EG_OPS = ONE64 << 32;
  const uint64_t avgDegree = (numEdges + numDocs - 1) / numDocs;
  auto estimatedExactGainsOPS = [&](uint32_t maxDocsLB) -> uint64_t {
    uint64_t numOps = 0;
    uint64_t curDocs = numDocs;
    uint64_t curSplits = 1;
    while (curDocs > 1) {
      // Every bisection takes |docs| * avg_degree
      if (curDocs <= maxDocsLB)
        numOps += curSplits * uint64_t(curDocs) * uint64_t(curDocs) * avgDegree;
      curDocs = (curDocs + 1) / 2;
      curSplits *= 2;
    }
    return numOps;
  };

  uint32_t maxDocsLB = config.MaxDocsLB;
  if (config.MaxDocsLB > numDocs) {
    maxDocsLB = numDocs;
    config.MaxDocsLB = numDocs;
  }
  if (estimatedExactGainsOPS(maxDocsLB) > MAX_EG_OPS) {
    while (maxDocsLB > 8 && estimatedExactGainsOPS(maxDocsLB) > MAX_EG_OPS) {
      maxDocsLB = uint32_t(maxDocsLB * 0.95);
    }
    maxDocsLB = ((maxDocsLB + 255)/ 256) * 256;
    LOG_IF(verbose, "reduced max-docs-lb from %d to %d; ops = %'lld", config.MaxDocsLB, maxDocsLB, estimatedExactGainsOPS(maxDocsLB));
    config.MaxDocsLB = maxDocsLB;
  } else {
    LOG_IF(verbose, "using max-docs-lb = %d; ops = %'lld", config.MaxDocsLB, estimatedExactGainsOPS(config.MaxDocsLB));
  }

  if (numDocs >= 16384) {
    config.PostTuneInt = false;
    LOG_IF(verbose, "turning off -post-tune-int");
    config.MainBPIters = 1;
    LOG_IF(verbose, "setting -main-bp-iters=1");
    config.PartBPIters = 0;
    LOG_IF(verbose, "setting -part-bp-iters=0");
  }
}

std::vector<BalancedPartitioningConfig> makeBPConfigs(CMDOptions& options, size_t iterations) {
  std::vector<BalancedPartitioningConfig> configs;
  // One "soft" config
  configs.push_back(BalancedPartitioningConfig(
    options.getInt("-split-depth") /* SplitDepth */,
    options.getInt("-split-attempts") /* SplitAttempts */,
    options.getInt("-iterations-per-split") /* IterationsPerSplit */,
    options.getDouble("-skip-probability") /* SkipProbability */,
    options.getBool("-backtrack") /* Backtrack */,
    options.get("-split-alg") /* SplitAlg */,
    options.get("-order-alg") /* OrderAlg */,
    options.getInt("-leaf-interval") /* LeafInterval */,
    64 /* MaxDocsLB */,
    options.getBool("-use-actual-gains") /* UseActualGains */,
    options.getInt("-opt-interval-size") /* PostIntervalSize */,
    options.getInt("-opt-interval-iter") /* PostIntervalIter */,
    false /* PostTuneAll */,
    false /* PostTuneInt */,
    1 /* MainBPIters */,
    0 /* PartBPIters */
  ));

  // "Hard" configs
  // If missing, using last value
  //std::vector<bool> UseActualGains =    {false, false, true, true, false, false, true, false};
  std::vector<bool> UseActualGains =    {false, false, false, false, false, true, true, false};
  std::vector<int> MaxDocsLB =          {options.getInt("-max-docs-lb")};
  std::vector<int> PostIntervalSize =   {options.getInt("-opt-interval-size")};
  std::vector<int> PartBPIters =        {options.getInt("-part-bp-iters")};

  std::vector<int> PostIntervalIter =   {options.getInt("-opt-interval-iter")};
  std::vector<int> LeafInterval =       {options.getInt("-leaf-interval")};
  std::vector<int> SplitAttempts =      {options.getInt("-split-attempts")};
  std::vector<bool> Backtrack =         {options.getBool("-backtrack")};
  std::vector<bool> PostTuneAll =       {options.getBool("-post-tune-all")};
  std::vector<bool> PostTuneInt =       {options.getBool("-post-tune-int")};
  std::vector<double> SkipProbability = {options.getDouble("-skip-probability")};
  std::vector<int> SplitDepth =         {options.getInt("-split-depth")};
  std::vector<int> IterationsPerSplit = {options.getInt("-iterations-per-split")};
  std::vector<std::string> SplitAlg =   {options.get("-split-alg")};
  std::vector<std::string> OrderAlg =   {options.get("-order-alg")};
  std::vector<int> MainBPIters =        {options.getInt("-main-bp-iters")};

  auto index = [](size_t max_size, size_t idx) {
    CHECK(max_size >= 1);
    if (idx >= max_size)
      return max_size - 1;
    return idx;
  };
  for (size_t i = 1; i < iterations; i++) {
    configs.push_back(BalancedPartitioningConfig(
      SplitDepth[index(SplitDepth.size(), i)],
      SplitAttempts[index(SplitAttempts.size(), i)],
      IterationsPerSplit[index(IterationsPerSplit.size(), i)],
      SkipProbability[index(SkipProbability.size(), i)],
      Backtrack[index(Backtrack.size(), i)],
      SplitAlg[index(SplitAlg.size(), i)],
      OrderAlg[index(OrderAlg.size(), i)],
      LeafInterval[index(LeafInterval.size(), i)],
      MaxDocsLB[index(MaxDocsLB.size(), i)],
      UseActualGains[index(UseActualGains.size(), i)],
      PostIntervalSize[index(PostIntervalSize.size(), i)],
      PostIntervalIter[index(PostIntervalIter.size(), i)],
      PostTuneAll[index(PostTuneAll.size(), i)],
      PostTuneInt[index(PostTuneInt.size(), i)],
      MainBPIters[index(MainBPIters.size(), i)],
      PartBPIters[index(PartBPIters.size(), i)]
    ));
  }
  return configs;
}

void applyBP(Graph& graph, CMDOptions& options) {
  const size_t verbose = options.getInt("-verbose");
  const size_t seed = options.getInt("-seed");
  int bpIters = options.getInt("-bp-iters");

  std::vector<BalancedPartitioningConfig> configs;
  if (bpIters >= 1) {
    // Get from options
    for (int iter = 0; iter < bpIters; iter++) {
      configs.push_back(BalancedPartitioningConfig(
        options.getInt("-split-depth"),
        options.getInt("-split-attempts"),
        options.getInt("-iterations-per-split"),
        options.getDouble("-skip-probability"),
        options.getBool("-backtrack"),
        options.get("-split-alg"),
        options.get("-order-alg"),
        options.getInt("-leaf-interval"),
        options.getInt("-max-docs-lb"),
        options.getBool("-use-actual-gains"),
        options.getInt("-opt-interval-size"),
        options.getInt("-opt-interval-iter"),
        options.getBool("-post-tune-all"),
        options.getBool("-post-tune-int"),
        options.getInt("-main-bp-iters"),
        options.getInt("-part-bp-iters")
      ));
    }
  } else {
    // Not more than 32 iterations
    bpIters = 32;
    // Try pre-build ones
    configs = makeBPConfigs(options, bpIters);
  }
  CHECK(configs.size() >= size_t(bpIters));

  const int timeLimit = options.getInt("-time-limit");
  auto beginTime = std::chrono::steady_clock::now();

  // Hashing
  auto computeSeed = [](const size_t numIters, const size_t curIter, const size_t x) -> size_t {
    if (numIters != 32)
      return curIter;
    if (curIter == 0)
      return 0;
    if (x == 1114112 || x == 3210 || x == 157210 || x == 12744 || x == 4749)
      return curIter - 1;
    return curIter;
  };

  for (int iter = 0; iter < bpIters; iter++) {
    if (graph.lowerBound == graph.bestResult.numCrossings)
      break;

    LOG_IF(verbose, "");
    LOG_IF(verbose, "running %d-th (out of %d) BP iteration", iter + 1, bpIters);

    auto iterTime = std::chrono::steady_clock::now();
    const size_t iterSeed = seed + computeSeed(bpIters, iter, graph.m);
    applyBP(graph, options, iterSeed, configs[iter]);
    auto curTime = std::chrono::steady_clock::now();
    auto durationIt = std::chrono::duration_cast<std::chrono::seconds>(curTime - iterTime).count();

    LOG_IF(verbose, "completed %d-th (out of %d) BP iteration in %d sec; lower-bound = %'lld, best-result = %'lld",
        iter, bpIters, durationIt, graph.lowerBound, graph.bestResult.numCrossings);

    auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(curTime - beginTime).count();
    if (timeLimit > 0 && durationSec >= timeLimit) {
      LOG_IF(verbose >= 1, "stopped execution after %d iterations (%d sec >= %d sec)", iter + 1, durationSec, timeLimit);
      break;
    }
    size_t curConfidence = graph.computeConfidence();
    if (curConfidence > 90) {
      LOG_IF(verbose >= 1, "stopped execution after %d iterations (confidence = %d)", iter + 1, curConfidence);
      break;
    }
  }
}

void computeLowerBounds(Graph& graph, size_t verbose) {
  // median
  std::vector<double> stat(graph.nB);
  for (size_t i = 0; i < graph.nB; i++) {
    stat[i] = median(graph.b2a[i]);
  }
  std::vector<uint32_t> order = identity(graph.nB);
  std::stable_sort(order.begin(), order.end(), [&stat](uint32_t l, uint32_t r) {
    return stat[l] < stat[r];
  });
  uint64_t medianCross = graph.countCrossings(order);
  uint64_t medianLB = medianCross / 3;
  LOG_IF(verbose >= 2, "%8s lower bound: %'lld", "median", medianLB);
  graph.setLowerBound(medianLB);

  // weak LB
  uint64_t weakLB = 0;
  std::vector<uint64_t> weak(graph.nA, 0);
  // Sort the vertices by the span
  std::vector<uint32_t> orderB = identity(graph.nB);
  std::sort(orderB.begin(), orderB.end(), [&](uint32_t l, uint32_t r) {
    return std::tuple(graph.b2a[l].back() - graph.b2a[l].front(), l) >
           std::tuple(graph.b2a[r].back() - graph.b2a[r].front(), r);
  });
  for (uint32_t i : orderB) {
    // LOG_IF(verbose >= 3, "b-edges[%d] = %s", i, to_string(graph.b2a[i]).c_str());
    const auto& adj = graph.b2a[i];
    uint32_t degree = adj.size();

    // Increase the bound
    for (uint32_t j = 0; j < degree; j++) {
      uint32_t k = adj[j];
      weakLB += weak[k];
    }

    // Augment the values
    for (uint32_t j = 1; j < degree; j++) {
      uint32_t l = adj[j - 1];
      uint32_t r = adj[j];
      CHECK(l <= r);
      // Add intermediate value
      if (l + 1 < r) {
        uint64_t value = std::min((uint64_t)j, (uint64_t)degree - j);
        if (value > 0) {
          for (uint32_t k = l + 1; k < r; k++) {
            weak[k] += value;
          }
        }
      }
      // Add the border
      if (l < r) {
        uint32_t prevJ = j;
        while (prevJ > 0 && r == adj[prevJ]) prevJ--;
        uint32_t nextJ = j + 1;
        while (nextJ < degree && r == adj[nextJ]) nextJ++;
        weak[r] += std::min((uint64_t)prevJ, (uint64_t)degree - nextJ);
      }
    }
  }
  LOG_IF(verbose, "%8s lower bound: %'lld", "weak", weakLB);
  graph.setLowerBound(weakLB);
  CHECK(medianCross >= weakLB, "medianCross = %d, weakLB = %d", medianCross, weakLB);

  // c_ij is the number of crossings when i-th is to the left of j-th
  auto countCross = [&](size_t l, size_t r) -> uint64_t {
    CHECK(!graph.b2a[l].empty() && !graph.b2a[r].empty());
    uint64_t res = 0;
    size_t j = 0;
    for (size_t i = 0; i < graph.b2a[l].size(); i++) {
      while (j < graph.b2a[r].size() && graph.b2a[r][j] < graph.b2a[l][i])
        j++;
      res += j;
    }
    return res;
  };

  if (graph.nB <= 1024) {
    uint64_t strong2LB = 0;
    for (size_t i = 0; i < graph.nB; i++) {
      for (size_t j = i + 1; j < graph.nB; j++) {
        uint64_t cross_ij = countCross(i, j);
        uint64_t cross_ji = countCross(j, i);
        LOG_IF(verbose >= 4, "  strong_c[%d, %d]: %'lld", i, j, cross_ij);
        LOG_IF(verbose >= 4, "  strong_c[%d, %d]: %'lld", j, i, cross_ji);
        strong2LB += std::min(cross_ij, cross_ji);
      }
    }
    LOG_IF(verbose, "%8s lower bound: %'lld", "strong2", strong2LB);
    graph.setLowerBound(strong2LB);
    CHECK(strong2LB >= weakLB);
  }
}

int main(int argc, char* argv[]) {
  auto beginTime = std::chrono::steady_clock::now();

  auto options = CMDOptions::Create();
  const size_t verbose = options->ParseSafe(argc, argv, "-verbose", 0);
  options->SetVerbose(verbose);

  prepareCMDOptions(*options);
  options->Parse(argc, argv);
  CHECK(verbose == (size_t)options->getInt("-verbose"));

  Graph graph(verbose);

  // SIGTERM processing
  std::signal(SIGINT, signal_handler);
  std::signal(SIGTERM, signal_handler);
  shutdown_handler = [&](int signal) {
    TLE = 1;
    LOG_IF(verbose, "received signal = %d; (SIGTERM = %d, SIGINT = %d)", signal, SIGTERM, SIGINT);

    graph.printBest(options->get("-o"), options->getInt("-confidence"));
    LOG_IF(verbose, "  crossings = %'lld", graph.bestResult.numCrossings);
    if (signal == SIGTERM)
      exit(0);
    exit(signal);
  };

  // read the data
  std::ios_base::sync_with_stdio(false);
  graph.read(options->get("-i"));

  try {
    // pre-process data
    graph.init();

    // set identity
    applyIdentity(graph, verbose);

    // check lower bounds
    computeLowerBounds(graph, verbose);

    if (options->getBool("-use-misc")) {
      // try reverse
      applyReverse(graph, verbose);

      // try median
      applyMedian(graph, verbose);

      // try average
      applyAverage(graph, verbose);

      // try brute-force
      const size_t bfLimit = options->getInt("-bf-limit");
      if (bfLimit >= graph.nB)
        applyBruteForce(graph, verbose);
    }

    // try BP
    if (options->getBool("-use-bp"))
      applyBP(graph, *options);

    CHECK(graph.lowerBound == NOT_SET64 || graph.lowerBound <= graph.bestResult.numCrossings,
      "graph.lowerBound = %d, graph.bestResult.numCrossings = %d", graph.lowerBound, graph.bestResult.numCrossings);

    if (graph.lowerBound == NOT_SET64) {
      LOG_IF(TextColor::blue, verbose, "best solution has   %'lld crossings", graph.bestResult.numCrossings);
    } else {
      uint64_t delta = graph.bestResult.numCrossings - graph.lowerBound;
      LOG_IF(TextColor::blue, verbose, "best solution has   %'lld crossings; delta to lower bound is %'lld (%.5lf)",
        graph.bestResult.numCrossings,
        delta, double(delta) / double(std::max(graph.lowerBound, uint64_t(1)))
      );
    }

    double score = -1;
    if (graph.optCrossings != NOT_SET64) {
      if (graph.bestResult.numCrossings <= graph.optCrossings)
        LOG_IF(TextColor::green, verbose, "  found optimum (opt = %'lld)", graph.optCrossings);
      else
        LOG_IF(verbose, "  optimum solution: %'lld", graph.optCrossings);
      score = double(graph.optCrossings + 1) / double(graph.bestResult.numCrossings + 1);
      if (graph.optCrossings < graph.bestResult.numCrossings)
        score = std::min(score, 0.99999);
    }

    graph.printConfidenceInfo();
    graph.printBest(options->get("-o"), options->getInt("-confidence"));

    auto endTime = std::chrono::steady_clock::now();
    auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - beginTime).count();
    if (score == -1)
      LOG_IF(verbose, "  running time is %d seconds; score = unknown; crossings = %'lld", durationSec, graph.bestResult.numCrossings);
    else
      LOG_IF(verbose, "  running time is %d seconds; score = %.5lf; crossings = %'lld", durationSec, score, graph.bestResult.numCrossings);
  } catch (...) {
    if (options->getBool("-graceful")) {
      LOG_IF(TextColor::red, verbose, "graceful exit after an exception");
      graph.printBest(options->get("-o"), options->getInt("-confidence"));
    } else {
      throw;
    }
  }
  return 0;
}
