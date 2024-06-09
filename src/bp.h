#pragma once

#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <unordered_set>

struct UtilitySignature;

/// An edge between a document and a utility vertex.
struct EdgeTy {
  uint32_t utility;
  uint32_t weight;

  explicit EdgeTy(uint32_t utility, uint32_t weight): utility(utility), weight(weight) {}
};

/// A document with edges to utility nodes. After merging duplicates, the object
/// may represent a group of duplicate (or very similar) documents whose ids
/// are stored in the corresponding field.
class Document {
public:
  Document(const Document &) = delete;
  Document(Document &&) = default;
  Document &operator=(const Document &) = delete;
  Document &operator=(Document &&) = default;

  explicit Document() {}

  struct IdOrder {
    bool operator()(const Document* L, const Document* R) const {
      return L->id() < R->id();
    }
  };

  struct BucketOrder {
    bool operator()(const Document* L, const Document* R) const {
      return std::make_pair(L->Bucket, L->id()) < std::make_pair(R->Bucket, R->id());
    }
  };

  struct SortIdxOrder {
    bool operator()(const Document* L, const Document* R) const {
      return L->SortIdx < R->SortIdx;
    }
  };

  struct MedianOrder {
    bool operator()(const Document* L, const Document* R) const {
      return std::tuple(L->AllMedian, L->AllAverage, L->id()) < std::tuple(R->AllMedian, R->AllAverage, R->id());
    }
  };

  struct LocalMedianOrder {
    bool operator()(const Document* L, const Document* R) const {
      return std::tuple(L->LocalMedian, L->LocalAverage, L->id()) < std::tuple(R->LocalMedian, R->LocalAverage, R->id());
    }
  };

  struct AverageOrder {
    bool operator()(const Document* L, const Document* R) const {
      return std::tuple(L->AllAverage, L->AllMedian, L->id()) < std::tuple(R->AllAverage, R->AllMedian, R->id());
    }
  };

  struct LocalAverageOrder {
    bool operator()(const Document* L, const Document* R) const {
      return std::tuple(L->LocalAverage, L->LocalMedian, L->id()) < std::tuple(R->LocalAverage, R->LocalMedian, R->id());
    }
  };

public:
  void init(uint32_t Id) { Ids.push_back(Id); }

  void init(const std::vector<uint32_t>& DocIds) {
    Ids.assign(DocIds.begin(), DocIds.end());
  }

  /// Get the (first) original document id.
  inline uint32_t id() const { return Ids.front(); }

  /// Get the list of original documents corresponding to the instance.
  const std::vector<uint32_t> ids() const { return Ids; }

  inline uint32_t degree() const { return Edges.size(); }

  uint64_t weightedDegree() const { return WeightedDegree; }

  const std::vector<EdgeTy>& getEdges() const { return Edges; }

  EdgeTy getEdge(uint32_t Index) const { return Edges[Index]; }

  void setEdges(const std::vector<EdgeTy>& UtilityNodes) {
    Edges.assign(UtilityNodes.begin(), UtilityNodes.end());
  }

  const std::vector<EdgeTy>& getOrigEdges() const { return OrigEdges; }

  void setOrigEdges(const std::vector<EdgeTy>& UtilityNodes) {
    Edges.assign(UtilityNodes.begin(), UtilityNodes.end());
    OrigEdges.assign(UtilityNodes.begin(), UtilityNodes.end());
    WeightedDegree = 0;
    for (const auto [_, W] : OrigEdges) {
      WeightedDegree += W;
    }
  }

  void restoreOrigEdges() { Edges = OrigEdges; }

  void incEdgeWeight(uint32_t Pos, uint32_t W) {
    Edges[Pos].weight += W;
    OrigEdges[Pos].weight += W;
    WeightedDegree += W;
  }

  uint64_t getSelfCrossings() const { return SelfCrossings; }

  void setSelfCrossings(uint64_t Crossings) { SelfCrossings = Crossings; }

  /// Original index of the document.
  uint32_t Index = uint32_t(-1);
  /// Local index of the document.
  uint32_t LocalIndex = uint32_t(-1);

  /// Document bucket assigned by balanced partitioning.
  uint32_t Bucket = uint32_t(-1);

  /// Move gains.
  int64_t MoveGain = 0;
  std::vector<uint32_t> AdjLB;
  uint32_t LastMovedIter = uint32_t(-1);

  /// SCCs.
  std::vector<uint32_t> SuccLB;
  std::vector<uint32_t> PredLB;

  /// Hash code of the document based on its content.
  double AllMedian = 0;
  double AllAverage = 0;
  double LocalMedian = 0;
  double LocalAverage = 0;
  /// Used for sorting.
  uint32_t SortIdx = 0;
  /// Used in tuning.
  uint32_t BestBucket = uint32_t(-1);

private:
  /// Document ids of all (duplicate) documents corresponding to the instance.
  std::vector<uint32_t> Ids;

  /// (current) Adjacent utility nodes of the document.
  std::vector<EdgeTy> Edges;

  /// (original) Adjacent utility nodes of the document.
  std::vector<EdgeTy> OrigEdges;

  /// The number of incident (non-merged) edges.
  uint64_t WeightedDegree = 0;

  /// The number of self-crossings.
  uint64_t SelfCrossings = 0;
};

using IterRange = std::pair<const std::vector<Document *>::iterator,
                            const std::vector<Document *>::iterator>;

enum class OrderAlgT { INPUT, RANDOM, MEDIAN, AVERAGE, LOCAL_MEDIAN, LOCAL_AVERAGE, NONE, MEDIAN_OPT, AVERAGE_OPT };
enum class SplitAlgT { INPUT, RANDOM, MEDIAN, AVERAGE };
enum class OptLevelT { O0, O1, O2, O3, O4 };

/// Algorithm parameters; default values are tuned on real-world instances.
struct BalancedPartitioningConfig {
  explicit BalancedPartitioningConfig(
    uint32_t SplitDepth,
    uint32_t SplitAttempts,
    uint32_t IterationsPerSplit,
    double SkipProbability,
    bool Backtrack,
    const std::string& SplitAlg,
    const std::string& OrderAlg,
    uint32_t LeafInterval,
    uint32_t MaxDocsLB,
    bool UseActualGains,
    uint32_t PostIntervalSize,
    uint32_t PostIntervalIter,
    bool PostTuneAll,
    bool PostTuneInt,
    uint32_t MainBPIters,
    uint32_t PartBPIters
  ):
    SplitDepth(SplitDepth),
    SplitAttempts(SplitAttempts),
    IterationsPerSplit(IterationsPerSplit),
    SkipProbability(SkipProbability),
    Backtrack(Backtrack),
    LeafInterval(LeafInterval),
    MaxDocsLB(MaxDocsLB),
    UseActualGains(UseActualGains),
    PostIntervalSize(PostIntervalSize),
    PostIntervalIter(PostIntervalIter),
    PostTuneAll(PostTuneAll),
    PostTuneInt(PostTuneInt),
    MainBPIters(MainBPIters),
    PartBPIters(PartBPIters)
  {
    setSplitAlg(SplitAlg);
    setOrderAlg(OrderAlg);
  }

  /// The depth of the recursive bisection.
  uint32_t SplitDepth = 18;

  /// The maximum attempts to perform a split.
  uint32_t SplitAttempts = 1;

  /// The maximum number of bp iterations per split.
  uint32_t IterationsPerSplit = 40;

  /// The probability for a vertex to skip a move from its current bucket to
  /// another bucket; it often helps to escape from a local optima.
  double SkipProbability = 0.01;

  /// Use backtracking.
  bool Backtrack = true;

  /// Split strategy.
  SplitAlgT SplitAlg = SplitAlgT::RANDOM;
  void setSplitAlg(const std::string& alg);

  /// Order strategy.
  OrderAlgT OrderAlg = OrderAlgT::LOCAL_MEDIAN;
  void setOrderAlg(const std::string& alg);

  /// For full ordering of leaves.
  uint32_t LeafInterval = 18;

  /// Maximum number of documents for lower bounds.
  uint32_t MaxDocsLB = 6144;

  /// Using up-to-date move gains for swapping.
  bool UseActualGains = false;

  /// The size of intervals for post-processing of reordered documents.
  uint32_t PostIntervalSize = 16;

  /// The size of intervals for post-processing of reordered documents.
  uint32_t PostIntervalIter = 15;

  /// Whether to apply fine-tuning at post-processing.
  bool PostTuneAll = true;

  /// Whether to apply fine-tuning at post-processing.
  bool PostTuneInt = true;

  /// The maximum number of simple-BP repetitions.
  uint32_t MainBPIters = 1;

  /// The maximum number of partial-BP repetitions.
  uint32_t PartBPIters = 0;
};

void adjustBPConfig(int verbose, BalancedPartitioningConfig& config, uint32_t numDocs, uint32_t numEdges);

/// Recursive balanced graph partitioning algorithm.
///
/// The algorithm is used to find an ordering of Documents while optimizing
/// a specified objective. The algorithm uses recursive bisection; it starts
/// with a collection of unordered documents and tries to split them into
/// two sets (buckets) of equal cardinality. Each bisection step is comprised of
/// iterations that greedily swap the documents between the two buckets while
/// there is an improvement of the objective. Once the process converges, the
/// problem is divided into two sub-problems of half the size, which are
/// recursively applied for the two buckets. The final ordering of the documents
/// is obtained by concatenating the two (recursively computed) orderings.
template <class LBDataType>
class BalancedPartitioning {
  using SignaturesType = std::vector<UtilitySignature>;

private:
  BalancedPartitioning(const BalancedPartitioning &) = delete;
  BalancedPartitioning &operator=(const BalancedPartitioning &) = delete;

public:
  explicit BalancedPartitioning(const BalancedPartitioningConfig& Config,
                                const bool UseDenseLBData);

  /// Initialize LB data before reordering.
  void initializeLB(std::vector<Document *>& Documents);

  /// Run recursive graph partitioning that optimizes a given objective.
  void run(std::vector<Document *>& Documents);

  /// Tune all intervals using a given optimization level.
  void tuneSolution(std::vector<Document *>& Documents,
                    std::function<void(const std::string& name)> progressFunc,
                    std::function<bool()> convergedFunc,
                    const OptLevelT OptLevel,
                    const size_t verbose);

  /// Tune all intervals using a given optimization level.
  void tuneLargeSolution(std::vector<Document *>& Documents,
                    std::function<void(const std::string& name)> progressFunc,
                    std::function<bool()> convergedFunc,
                    const OptLevelT OptLevel,
                    const size_t verbose);

  /// Compute a lower bounds based on opt-intervals.
  void computeLowerBound(std::vector<Document *>& Documents,
                         uint32_t IntervalSize);

  /// Return the lower bound based on opt-intervals.
  uint64_t getLowerBound() const { return LowerBound; };

  void setVerbose(size_t v) { verbose = v; }

  size_t getVerbose() const { return verbose; }

  void setSeed(size_t s) { seed = s; }

  size_t getSeed() const { return seed; }

 private:
  /// Initialize documents before reordering.
  void initialize(std::vector<Document *>& Documents);

  /// Run a recursive bisection of a given list of documents, where
  ///  - 'RecDepth' is the current depth of recursion
  ///  - 'RootBucket' is the initial bucket of the documents
  ///  - the assigned buckets are the range [Offset, Offset + NumDocuments)
  void bisect(const std::vector<Document *>::iterator &DocumentBegin,
              const std::vector<Document *>::iterator &DocumentEnd,
              const uint32_t RecDepth, const uint32_t RootBucket,
              const uint32_t Offset, const bool SplitBySCCs);

  /// Run bisection iterations.
  /// Returns true iff a progress has been made.
  void runIterations(SignaturesType &Signatures,
                     const std::vector<Document *>::iterator &DocumentBegin,
                     const std::vector<Document *>::iterator &DocumentEnd,
                     const uint32_t RecDepth);

  /// Run a bisection iteration to improve the optimization goal.
  /// Returns the total number of moved documents.
  uint32_t runIteration(SignaturesType &Signatures,
                        const std::vector<Document *>::iterator &DocumentBegin,
                        const std::vector<Document *>::iterator &DocumentEnd,
                        const uint32_t Iter);

  /// Try to move a document from one bucket to another.
  /// Return true iff the document is moved.
  bool moveDocument(Document *Doc, SignaturesType &Signatures) const;

  void restoreOrigUtilities(const std::vector<Document*>::iterator& DocumentBegin,
                            const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Initialize utility node signatures.
  void initializeSignatures(SignaturesType &Signatures,
                            const std::vector<Document *>::iterator &DocumentBegin,
                            const std::vector<Document *>::iterator &DocumentEnd) const;

  /// Initialize utility node signatures for gains.
  void initializeSignatureFastGains(SignaturesType &Signatures,
                                    const std::vector<Document *>::iterator &DocumentBegin,
                                    const std::vector<Document *>::iterator &DocumentEnd) const;

  /// Update utility node signatures.
  void computeFastMoveGains(SignaturesType &Signatures,
                            const std::vector<Document *>::iterator &DocumentBegin,
                            const std::vector<Document *>::iterator &DocumentEnd) const;

  /// (Slow) move gains.
  void computeStrongMoveGains(const std::vector<Document *>::iterator &DocumentBegin,
                              const std::vector<Document *>::iterator &DocumentEnd) const;

  /// (Slow) move gain for a specific document.
  void updateStrongMoveGains(Document *Doc, uint32_t DocIndex,
                             const std::vector<Document *>::iterator &DocumentBegin,
                             const std::vector<Document *>::iterator &DocumentEnd) const;

  /// Strong move gain
  inline int64_t strongMoveGain(const uint64_t LB1, const uint64_t LB2,
                                const uint32_t Doc1Bucket, const uint32_t Doc2Bucket) const;

  /// Split all the documents into 2 buckets, StartBucket and StartBucket + 1.
  /// The method is used for an initial assignment before a bisection step.
  void split(const std::vector<Document *>::iterator &DocumentBegin,
             const std::vector<Document *>::iterator &DocumentEnd,
             const uint32_t StartBucket) const;

  /// Order the list of documents by assigning buckets in the range
  /// [StartBucket + documentBegin, StartBucket + documentEnd).
  /// The method is used for assigning buckets when the number of
  /// documents is small (to truncate the bisection tree).
  void order(const std::vector<Document *>::iterator &DocumentBegin,
             const std::vector<Document *>::iterator &DocumentEnd,
             uint32_t StartBucket, OrderAlgT OrderAlgorithm) const;

  std::pair<uint64_t, OrderAlgT> orderForBacktrack(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    const uint32_t StartBucket, const uint32_t NumUtilities,
    std::vector<uint32_t>& OrgIdOrder) const;

  /// Count crossings for the documents.
  uint64_t countCrossings(
      const std::vector<Document*>::iterator& DocumentBegin,
      const std::vector<Document*>::iterator& DocumentEnd,
      const uint32_t NumUtilities) const;

  /// Count crossings for the documents.
  uint64_t countCrossings(
      const std::vector<Document*>::iterator& DocumentBegin,
      const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Compute an interval-based lower bound.
  uint64_t computeLowerBound(
      const std::vector<Document*>::iterator& DocumentBegin,
      const std::vector<Document*>::iterator& DocumentEnd,
      const uint32_t IntervalSize) const;

  /// Estimate the objective value of partitioning (e.g., the expected number of crossings).
  uint64_t objective(SignaturesType &Signatures,
                     const std::vector<Document*>::iterator& DocumentBegin,
                     const std::vector<Document*>::iterator& DocumentEnd) const;

  void collectBestBuckets(const std::vector<Document*>::iterator& DocumentBegin,
                          const std::vector<Document*>::iterator& DocumentEnd) const;

  void assignBestBuckets(const std::vector<Document*>::iterator& DocumentBegin,
                         const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Compute a pair (current crossings, optimal crossings) for a given interval.
  std::pair<uint64_t, uint64_t> computeCurOptCrossings(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Check if the order of documents is (provably) optimal.
  bool isOptimal(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Sort the documents with respect to a given permutation.
  void sortByOrder(
      const std::vector<Document*>::iterator& DocumentBegin,
      const std::vector<Document*>::iterator& DocumentEnd,
      const std::vector<uint32_t>& Order) const;

  /// Split and reorder the range of documents by their SCCs.
  uint64_t splitAndReorderBySCCs(
      const std::vector<Document *>::iterator &DocumentBegin,
      const std::vector<Document *>::iterator &DocumentEnd,
      std::vector<IterRange>& Result) const;

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// LBs ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  void initializeLBMatrix(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const;

  void initLB(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const;

  bool isInitializedLB(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const;

  //////////////////////////////////////////////////////////////////////////////
  /////////////////////////// Fine-Tuning //////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /// Tune all intervals using a given optimization level.
  uint64_t tuneAllIntervals(
    std::vector<Document *>& Documents,
    const std::vector<std::pair<uint32_t, uint32_t>>& Intervals,
    std::function<void(const std::string& name)> progressFunc,
    const OptLevelT OptLevel,
    const size_t verbose);

  /// Tune a given interval using a given optimization level.
  uint64_t tuneInterval(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool SplitBySCCs, const OptLevelT Level,
    const size_t verbose) const;

  /// Sort the documents optimally via DP. Returns true if there is progress.
  uint64_t sortIntervalOpt(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    bool IsIncremental) const;

  /// Post-processing sorting intervals.
  uint64_t sortIntervalsOpt(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t IntervalSize,
    const uint32_t IntervalIter,
    const OptLevelT OptLevel,
    const size_t verbose) const;

  /// Merge two intervals optimally. Returns true if there is progress.
  uint64_t mergeTwoIntervals(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<bool(uint32_t, const Document*)> IsFirstFilter) const;

  /// Merge two intervals optimally. Returns true if there is progress.
  uint64_t mergeTwoIntervals(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentMid,
    const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Merge two degree-based intervals optimally.
  uint64_t mergeTwoIntervalsByDegree(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    uint32_t MaxDegree) const;

  /// Merge two odd/even-based intervals optimally.
  uint64_t mergeTwoIntervalsOddEven(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const;

  /// Merge interval pairs optimally.
  uint64_t mergeIntervalPairs(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t MaxDiffSize) const;

  /// Merge two intervals induced by the given filter optimally.
  uint64_t adjustWithFilter(
      std::function<bool(uint32_t, const Document*)> Filter,
      const std::string FilterName) const;

  /// Place all singletons optimally.
  uint64_t adjustSingle(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t FirstIndex,
    std::vector<uint64_t>& SumL,
    std::vector<uint64_t>& SumR) const;

  /// Place all singletons optimally.
  uint64_t adjustSingle(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const OptLevelT OptLevel) const;

  /// Random search for a given set of documents.
  uint64_t randomSearchIncremental(
    std::vector<Document *>& Documents,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t IntervalSize,
    const uint32_t IntervalIter) const;

  /// Random search for a given set of documents.
  uint64_t randomSearchIncremental(
      const std::vector<Document*>::iterator& DocumentBegin,
      const std::vector<Document*>::iterator& DocumentEnd,
      const uint32_t MaxAttempts,
      const bool SplitBySCCs) const;

  /// Random search for a given set of documents.
  uint64_t randomSearchWithSwaps(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t NumSwaps,
    const uint32_t MaxIters,
    const uint32_t MaxItersWithoutProgress) const;

  /// Random search for a given set of documents.
  uint64_t randomSearch(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool SplitBySCCs,
    std::function<void(const std::vector<Document*>::iterator&, const std::vector<Document*>::iterator&, const uint32_t)> genOrderFunc,
    std::function<bool(const std::vector<uint64_t>&)> convergeFunc) const;

 private:
  /// Algorithm parameters.
  const BalancedPartitioningConfig &Config;

  /// Input documents that shall be reordered by the algorithm.
  uint32_t MaxUtility = 0;

  /// For logging
  size_t verbose = 0;
  /// Random seed
  size_t seed = 0;
  /// Random-number generator.
  mutable std::mt19937 RNG;

  /// Pre-allocated segment-tree.
  mutable std::vector<uint32_t> stree;

  /// Pre-allocated memory for sorting intervals.
  static constexpr uint32_t MaxDocsDP = 18;
  static constexpr uint32_t MaxDPSize = uint32_t(1) << MaxDocsDP;
  mutable std::vector<uint32_t> DP;
  mutable std::vector<uint16_t> DPLast;
  mutable std::vector<uint16_t> DPFirst;
  mutable std::vector<std::vector<uint16_t>> DPPrevs;
  mutable LBDataType LB_rev[MaxDocsDP][MaxDocsDP];
  mutable std::vector<std::vector<uint32_t>> DPMasks;

  /// Pre-allocated memory for interval merging
  mutable std::vector<std::vector<uint64_t>> MergeDP;
  mutable std::vector<std::vector<uint64_t>> LBSum1;
  mutable std::vector<std::vector<uint64_t>> LBSum2;

  /// For lower bounds
  uint64_t LowerBound = uint64_t(-1);

  /// Current partitioning
  bool UseStrongMoveGains = false;
  uint32_t LeftBucket = uint32_t(-1);
  uint32_t RightBucket = uint32_t(-1);

  /// Sparse/Dense LB data
  const bool UseDenseLBData;
  /// LB matrix
  mutable std::vector<std::vector<LBDataType>> LBData;
  /// Cached doc indices (for DenseLB).
  mutable std::unordered_set<uint32_t> CachedDocs;
};

/// Signature of a Utility Node utilized in a bisection step, that is, the
/// number of incident documents in the two buckets.
struct UtilitySignature {
  UtilitySignature(const UtilitySignature&) = delete;
  UtilitySignature(UtilitySignature&&) = default;
  UtilitySignature& operator=(const UtilitySignature&) = delete;
  UtilitySignature& operator=(UtilitySignature&&) = default;

  explicit UtilitySignature(uint32_t LeftCount = 0, uint32_t RightCount = 0)
      : LeftCount(LeftCount), RightCount(RightCount) {}

 public:
  /// The number of documents in the left bucket.
  uint32_t LeftCount;
  /// The number of documents in the right bucket.
  uint32_t RightCount;

  /// Cached cost of moving a document from left to right bucket.
  double CachedCostLR = 0;
  /// Aux values for weak gains.
  int64_t YRB = 0;
  int64_t YLB = 0;
  int64_t XRB = 0;
  int64_t XLB = 0;
};

/// Update document adjacency lists. Returns the number of utilities adjacent
/// to the given set of documents (that is, MaxUtilitity + 1)
uint32_t updateUtilities(const std::vector<Document*>::iterator& DocumentBegin,
                         const std::vector<Document*>::iterator& DocumentEnd,
                         bool UpdateOrigEdges=false);

/// Get the sorted document ids.
std::vector<uint32_t> extractOrderedIds(const std::vector<Document *>& Docs);

/// Check if buckets are correctly ordered.
bool bucketsCorrect(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool StrictOrder=false);

/// Check if the optimization coverged.
bool attemptsConverged(
    const std::vector<uint64_t>& IterObjective,
    const uint32_t kMaxIters,
    const uint32_t kMaxItersWithoutProgress);

/// Collect current buckets to a specified vector.
void collectBuckets(const std::vector<Document *>& Documents,
                    std::vector<std::pair<uint32_t, uint32_t>>& Buckets);

/// Assign buckets from a vector to current buckets.
void assignBuckets(std::vector<Document *>& Documents,
                   const std::vector<std::pair<uint32_t, uint32_t>>& Buckets);
