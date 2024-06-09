#include "bp.h"
#include "logging.h"
#include "common.h"

bool attemptsConverged(
    const std::vector<uint64_t>& IterObjective,
    const uint32_t kMaxIters,
    const uint32_t kMaxItersWithoutProgress) {
  if (IterObjective.size() > kMaxIters)
    return true;
  CHECK(IterObjective.size() >= 1);
  uint64_t BestCrossings = IterObjective[0];
  uint32_t ItersWithoutProgress = 0;
  for (uint32_t Iter = 1; Iter < IterObjective.size(); Iter++) {
    const uint64_t NewCrossings = IterObjective[Iter];
    if (NewCrossings < BestCrossings) {
      BestCrossings = NewCrossings;
      ItersWithoutProgress = 0;
    } else {
      ItersWithoutProgress++;
    }
    if (ItersWithoutProgress >= kMaxItersWithoutProgress)
      return true;
  }
  uint32_t kMaxEqualIters = 5;
  if (IterObjective.size() >= kMaxEqualIters) {
    uint32_t NumBestAttempts = 0;
    for (uint32_t I = 0; I < kMaxEqualIters; I++) {
      if (IterObjective[IterObjective.size() - I - 1] == BestCrossings) {
        NumBestAttempts++;
      }
    }
    if (NumBestAttempts == kMaxEqualIters)
      return true;
  }
  return false;
}

std::vector<std::pair<uint32_t, uint32_t>>
splitIntoIntervals(uint32_t NumDocuments, uint32_t MaxDocsLB, uint32_t Ratio) {
  CHECK(MaxDocsLB > 0);
  if (NumDocuments <= MaxDocsLB) {
    return {{0, NumDocuments}};
  }
  CHECK(Ratio > 1);

  std::vector<std::pair<uint32_t, uint32_t>> Res;
  uint32_t Delta = std::max(uint32_t(1), MaxDocsLB / Ratio);
  //uint32_t Delta = std::max(uint32_t(1), MaxDocsLB / 4);
  for (uint32_t I = 0; I < NumDocuments; I += Delta) {
    uint32_t EI = std::min(NumDocuments, I + MaxDocsLB);
    Res.push_back({I, EI});
    if (EI == NumDocuments)
      break;
  }

  return Res;
}

/// O0 == tuneAllIntervals
/// O1 == O0 + randomSearchIncremental
/// O2 == O1 + pre-randomSearchWithSwaps
/// O3 == O1 + pre-randomSearchWithSwaps + post-randomSearchWithSwaps
template <class T>
void BalancedPartitioning<T>::tuneSolution(
    std::vector<Document *>& Documents,
    std::function<void(const std::string& name)> progressFunc,
    std::function<bool()> convergedFunc,
    const OptLevelT OptLevel,
    const size_t verbose) {
  if (convergedFunc()) return;

  // For determinism
  RNG.seed(seed);

  std::vector<std::pair<uint32_t, uint32_t>> Intervals = splitIntoIntervals(Documents.size(), Config.MaxDocsLB, /* Ratio */ 2);
  tuneAllIntervals(Documents, Intervals, progressFunc, OptLevel, verbose);
  if (convergedFunc()) return;
  if (OptLevel == OptLevelT::O0)
    return;

  uint64_t SavedCrossingsBySwaps = 0;
  if (Config.PostTuneInt && OptLevel >= OptLevelT::O1) {
    const uint32_t numDocs = Documents.size();
    const std::vector<Document*>::iterator DocBegin = Documents.begin();
    SavedCrossingsBySwaps += randomSearchWithSwaps(DocBegin, DocBegin + numDocs / 3,                   progressFunc, 3, /* MaxIters */ 30, /* MaxItersWithoutProgress */ 15);
    SavedCrossingsBySwaps += randomSearchWithSwaps(DocBegin + numDocs / 3, DocBegin + 2 * numDocs / 3, progressFunc, 3, /* MaxIters */ 30, /* MaxItersWithoutProgress */ 15);
    SavedCrossingsBySwaps += randomSearchWithSwaps(DocBegin + 2 * numDocs / 3, DocBegin + numDocs,     progressFunc, 3, /* MaxIters */ 30, /* MaxItersWithoutProgress */ 15);
  }

  if (Config.PostTuneInt && OptLevel >= OptLevelT::O2 /*&& SavedCrossingsBySwaps > 0*/) {
    SavedCrossingsBySwaps += randomSearchWithSwaps(Documents.begin(), Documents.end(), progressFunc, 3, /* MaxIters */ 75, /* MaxItersWithoutProgress */ 20);
    if (SavedCrossingsBySwaps > 0)
    SavedCrossingsBySwaps += randomSearchWithSwaps(Documents.begin(), Documents.end(), progressFunc, 9, /* MaxIters */ 75, /* MaxItersWithoutProgress */ 10);
    if (convergedFunc()) return;
  }

  uint64_t SavedCrossingsByInt = 0;
  if (Config.PostTuneInt && OptLevel >= OptLevelT::O1) {
    std::vector<uint32_t> IntervalSizes = {50, 75, 100, 150};
    for (size_t I = 0; I < IntervalSizes.size(); I++) {
      SavedCrossingsByInt += randomSearchIncremental(Documents, progressFunc, IntervalSizes[I], /* MaxIters */ 100);
      if (convergedFunc()) return;
    }
  }

  if (Config.PostTuneInt && OptLevel >= OptLevelT::O3 && (SavedCrossingsByInt > 0 || SavedCrossingsBySwaps > 0)) {
    randomSearchWithSwaps(Documents.begin(), Documents.end(), progressFunc, 3, /* MaxIters */ 25, /* MaxItersWithoutProgress */ 20);
    randomSearchWithSwaps(Documents.begin(), Documents.end(), progressFunc, 9, /* MaxIters */ 25, /* MaxItersWithoutProgress */ 10);
    if (convergedFunc()) return;
  }
}

template <class T>
void BalancedPartitioning<T>::tuneLargeSolution(
    std::vector<Document *>& Documents,
    std::function<void(const std::string& name)> progressFunc,
    std::function<bool()> convergedFunc,
    const OptLevelT OptLevel,
    const size_t verbose) {
  if (convergedFunc()) return;

  // For determinism
  RNG.seed(seed);

  std::vector<std::pair<uint32_t, uint32_t>> Intervals;

  Intervals = splitIntoIntervals(Documents.size(), Config.MaxDocsLB, /* Ratio */ 3);
  tuneAllIntervals(Documents, Intervals, progressFunc, OptLevel, verbose);
  if (convergedFunc()) return;
  progressFunc("bp_" + std::to_string(seed) + "_tunelarge");

  Intervals = splitIntoIntervals(Documents.size(), Config.MaxDocsLB, /* Ratio */ 2);
  tuneAllIntervals(Documents, Intervals, progressFunc, OptLevel, verbose);
}

/// O0 == adjustSingle only
/// O1 == IntervalSize=11, IntervalIters=1
/// O2 == IntervalSize=11, IntervalIters=15
/// O3 == IntervalSize=14, IntervalIters=5
/// O4 == IntervalSize=16, IntervalIters=15
template <class T>
uint64_t BalancedPartitioning<T>::tuneInterval(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool SplitBySCCs,
    const OptLevelT OptLevel,
    const size_t verbose) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);

  CHECK(0 < N && N <= Config.MaxDocsLB);
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));

  // // Stop early if already optimal
  // if (isOptimal(DocumentBegin, DocumentEnd))
  //   return 0;

  // Stop early for a few documents
  if (N <= Config.LeafInterval) {
    uint64_t SavedCrossings = sortIntervalOpt(DocumentBegin, DocumentEnd, false);
    return SavedCrossings;
  }

  // Split by SCCs
  if (SplitBySCCs) {
    std::vector<IterRange> Ranges;
    uint64_t SavedCrossings = splitAndReorderBySCCs(DocumentBegin, DocumentEnd, Ranges);
    for (const auto [DocBegin, DocEnd] : Ranges) {
      if (std::distance(DocBegin, DocEnd) > 1) {
        SavedCrossings += tuneInterval(DocBegin, DocEnd, false, OptLevel, verbose);
      }
    }
    return SavedCrossings;
  }

  // Level-O0 optimizations
  uint64_t SavedCrossingsPre = 0;
  SavedCrossingsPre += mergeTwoIntervalsByDegree(DocumentBegin, DocumentEnd, 1);
  SavedCrossingsPre += mergeTwoIntervalsOddEven(DocumentBegin, DocumentEnd);
  SavedCrossingsPre += adjustSingle(DocumentBegin, DocumentEnd, OptLevel);

  if (OptLevel == OptLevelT::O0) {
    return SavedCrossingsPre;
  }

  CHECK(OptLevel >= OptLevelT::O1);
  CHECK(Config.PostIntervalSize >= 4);
  const uint32_t IntervalSize =
      OptLevel == OptLevelT::O1 ? std::min(uint32_t(11), Config.PostIntervalSize - 5) :
      OptLevel == OptLevelT::O2 ? std::min(uint32_t(11), Config.PostIntervalSize - 5) :
      OptLevel == OptLevelT::O3 ? std::min(uint32_t(14), Config.PostIntervalSize - 3) :
      Config.PostIntervalSize;
  const uint32_t IntervalIter =
      OptLevel == OptLevelT::O1 ? 1 :
      OptLevel == OptLevelT::O2 ? 15 :
      OptLevel == OptLevelT::O3 ? 5 :
      Config.PostIntervalIter;

  // Sorting intervals
  auto progressFunc = [](const std::string& name) {};
  uint64_t SavedCrossingsSort = sortIntervalsOpt(DocumentBegin, DocumentEnd, progressFunc, IntervalSize, IntervalIter, OptLevel, verbose);

  if (SavedCrossingsSort == 0) {
    return SavedCrossingsPre;
  }

  // Repeat this just in case
  uint64_t SavedCrossingsPost = 0;
  SavedCrossingsPost += mergeTwoIntervalsByDegree(DocumentBegin, DocumentEnd, 1);
  SavedCrossingsPost += mergeTwoIntervalsOddEven(DocumentBegin, DocumentEnd);
  SavedCrossingsPost += adjustSingle(DocumentBegin, DocumentEnd, OptLevel);
  if (OptLevel >= OptLevelT::O2) {
    SavedCrossingsPost += mergeIntervalPairs(DocumentBegin, DocumentEnd, 150);
  }

  if (SavedCrossingsPost == 0 || OptLevel == OptLevelT::O1) {
    return SavedCrossingsPre + SavedCrossingsSort + SavedCrossingsPost;
  }

  CHECK(OptLevel >= OptLevelT::O2);
  if (SavedCrossingsPost > 0) {
    SavedCrossingsSort += sortIntervalsOpt(DocumentBegin, DocumentEnd, progressFunc, IntervalSize, IntervalIter, OptLevel, verbose);
  }

  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));
  return SavedCrossingsPre + SavedCrossingsSort + SavedCrossingsPost;
}

template <class T>
uint64_t BalancedPartitioning<T>::tuneAllIntervals(
    std::vector<Document *>& Documents,
    const std::vector<std::pair<uint32_t, uint32_t>>& Intervals,
    std::function<void(const std::string& name)> progressFunc,
    const OptLevelT OptLevel,
    const size_t verbose) {
  using namespace std::chrono;

  CHECK(Config.MaxDocsLB > 0);
  auto beginTime = steady_clock::now();

  LOG_IF(verbose, "tune up full BP solution; got %d interval(s)", Intervals.size());

  uint64_t SavedCrossings = 0;
  uint32_t ImprovedIntervals = 0;
  for (const auto& [BI, EI] : Intervals) {
    const std::vector<Document*>::iterator DocumentBegin = Documents.begin() + BI;
    const std::vector<Document*>::iterator DocumentEnd = Documents.begin() + EI;

    // Pre-compute LBs
    initLB(DocumentBegin, DocumentEnd);
    // Need to optimize anything?
    if (isOptimal(DocumentBegin, DocumentEnd))
      continue;

    const size_t verboseInt = Intervals.size() == 1 ? verbose : 0;
    SavedCrossings += tuneInterval(DocumentBegin, DocumentEnd, true, OptLevel, verboseInt);
    if (SavedCrossings > 0) {
      ImprovedIntervals++;
      // Saving progress every 10 intervals
      if ((ImprovedIntervals + 1) % 10 == 0)
        progressFunc("bp_" + std::to_string(seed) + "_tuneall_" + std::to_string(ImprovedIntervals));
    }
  }
  LOG_IF(verbose >= 1, "  improved %d intervals out of %d", ImprovedIntervals, Intervals.size());

  // Re-compute the lower bound
  if (OptLevel >= OptLevelT::O2 && Documents.size() <= Config.MaxDocsLB) {
    computeLowerBound(Documents, Config.PostIntervalSize);
  }

  auto endTime = steady_clock::now();
  auto durationSec = duration_cast<seconds>(endTime - beginTime).count();
  LOG_IF(verbose >= 1, "  removed %'lld crossings by full tune up; running time is %d seconds", SavedCrossings, durationSec);

  if (OptLevel >= OptLevelT::O2) {
    // Always apply the check, as it can improve the lower bound
    progressFunc("bp_" + std::to_string(seed) + "_tuneall");
  }

  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::randomSearch(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool SplitBySCCs,
    std::function<void(const std::vector<Document*>::iterator&, const std::vector<Document*>::iterator&, const uint32_t)> genOrderFunc,
    std::function<bool(const std::vector<uint64_t>&)> convergeFunc) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);

  CHECK(0 < N && N <= Config.MaxDocsLB);
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  // Do not run for small intervals
  if (N <= Config.LeafInterval) {
    uint64_t SavedCrossings = sortIntervalOpt(DocumentBegin, DocumentEnd, false);
    return SavedCrossings;
  }

  // Split by individual components
  if (SplitBySCCs) {
    std::vector<IterRange> Ranges;
    uint64_t SavedCrossings = splitAndReorderBySCCs(DocumentBegin, DocumentEnd, Ranges);
    CHECK(Ranges.size() >= 1);
    for (const auto [DocBegin, DocEnd] : Ranges) {
      CHECK(DocBegin != DocEnd);
      if (std::distance(DocBegin, DocEnd) > 1) {
        SavedCrossings += randomSearch(DocBegin, DocEnd, false, genOrderFunc, convergeFunc);
      }
    }
    return SavedCrossings;
  }

  uint64_t SavedCrossings = 0;

  // Check the lower bound
  const auto [CurCrossings, OptCrossings] = computeCurOptCrossings(DocumentBegin, DocumentEnd);
  // Stop if already optimal
  if (CurCrossings == OptCrossings)
    return SavedCrossings;
  // need this?
  // SavedCrossings = tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O2, 0);
  // const uint64_t PostIntervalSize = 11;
  // const uint64_t LowerBoundInt = computeLowerBound(DocumentBegin, DocumentEnd, PostIntervalSize);
  // if (CurCrossings <= LowerBoundInt) {
  //   LOG_IF(verbose >= 1, "  achieved lower-bound-int: %'lld; NumDocuments = %d", LowerBoundInt, N);
  //   return SavedCrossings;
  // }

  uint64_t BestCrossings = CurCrossings;
  collectBestBuckets(DocumentBegin, DocumentEnd);

  // first element is always the current objective
  std::vector<uint64_t> IterObjectives;
  IterObjectives.push_back(CurCrossings);
  uint32_t AppliedIters = 0;
  for (uint32_t Iter = 0; ; Iter++) {
    AppliedIters = Iter + 1;
    // Shuffle the documents
    assignBestBuckets(DocumentBegin, DocumentEnd);
    genOrderFunc(DocumentBegin, DocumentEnd, Iter);

    const auto [StartCrossings, _] = computeCurOptCrossings(DocumentBegin, DocumentEnd);
    // Re-order by adjustments
    const uint64_t Delta = tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O1, 0);
    CHECK(Delta <= StartCrossings);
    const uint64_t NewCrossings = StartCrossings - Delta;

    LOG_IF(verbose >= 3, "    NewCrossings = %4d at iteration %d with %d docs", NewCrossings, Iter, N);

    if (NewCrossings < BestCrossings) {
      LOG_IF(verbose >= 3, "    removed %4d crossings by random search at iteration %d with %d docs", BestCrossings - NewCrossings, Iter, N);
      collectBestBuckets(DocumentBegin, DocumentEnd);
      BestCrossings = NewCrossings;
      if (BestCrossings == OptCrossings)
        break;
    }

    // Check if the process converged
    IterObjectives.push_back(NewCrossings);
    if (convergeFunc(IterObjectives))
      break;
  }

  // Sort documents by the resulting (or original) best buckets
  assignBestBuckets(DocumentBegin, DocumentEnd);
  std::sort(DocumentBegin, DocumentEnd, Document::BucketOrder());
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));

  CHECK(BestCrossings <= CurCrossings);
  SavedCrossings += CurCrossings - BestCrossings;

  LOG_IF(verbose >= 2, "  removed %d (%'lld -> %'lld) crossings by tuning intervals for %d docs using %d iterations",
    SavedCrossings, CurCrossings, BestCrossings, N, AppliedIters
  );

  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::randomSearchIncremental(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t MaxAttempts,
    const bool SplitBySCCs) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(0 < N && N <= Config.MaxDocsLB);
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  // Do not run for small intervals
  if (N <= Config.LeafInterval) {
    uint64_t SavedCrossings = sortIntervalOpt(DocumentBegin, DocumentEnd, false);
    return SavedCrossings;
  }

  if (SplitBySCCs) {
    std::vector<IterRange> Ranges;
    uint64_t SavedCrossings = splitAndReorderBySCCs(DocumentBegin, DocumentEnd, Ranges);
    CHECK(Ranges.size() >= 1);
    for (const auto [DocBegin, DocEnd] : Ranges) {
      CHECK(DocBegin != DocEnd);
      if (std::distance(DocBegin, DocEnd) > 1) {
        SavedCrossings += randomSearchIncremental(DocBegin, DocEnd, MaxAttempts, false);
      }
    }
    return SavedCrossings;
  }

  // Apply swaps
  const uint32_t MinBucket = DocumentBegin[0]->Bucket;
  auto genOrderFunc = [&](
      const std::vector<Document*>::iterator& DocBegin,
      const std::vector<Document*>::iterator& DocEnd,
      const uint32_t Iter) {
    if (Iter == 0) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::MEDIAN);
    } else if (Iter == 1) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::AVERAGE);
    } else if (Iter == 2) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::INPUT);
    } else if (Iter <= 8) {
      // RANDOM_SWAPS
      const uint32_t NumSwaps = 4;
      assignBestBuckets(DocBegin, DocEnd);
      for (uint32_t J = 0; J < NumSwaps; J++) {
        uint32_t D1 = Rand::next(N, RNG);
        uint32_t D2 = Rand::next(N, RNG);
        if (D1 == D2)
          continue;
        std::swap(DocBegin[D1]->Bucket, DocBegin[D2]->Bucket);
      }
      std::sort(DocBegin, DocEnd, Document::BucketOrder());
    } else {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::RANDOM);
    }
  };
  // Convergence
  const uint32_t MaxItersWithoutProgress = 20;
  auto convergeFunc = [&](const std::vector<uint64_t>& IterObjective) {
    return attemptsConverged(
      IterObjective,
      MaxAttempts,
      MaxItersWithoutProgress
    );
  };

  const uint64_t SavedCrossings = randomSearch(DocumentBegin, DocumentEnd, false, genOrderFunc, convergeFunc);
  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::randomSearchIncremental(
    std::vector<Document *>& Documents,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t IntervalSize,
    const uint32_t RandomAttempts) const {
  const uint32_t N = Documents.size();
  if (IntervalSize > N || IntervalSize > Config.MaxDocsLB)
    return 0;

  auto beginTime = std::chrono::steady_clock::now();
  LOG_IF(verbose, "random search for intervals of size = %d using %d attempts", IntervalSize, RandomAttempts);
  uint64_t SavedCrossings = 0;

  // Build the intervals
  const uint32_t IntervalIters = 5;
  const uint32_t StepSize = std::max(IntervalSize / 4, uint32_t(15));

  std::vector<bool> NeedUpdate(N, true);
  std::vector<bool> NeedUpdateNext(N, false);
  for (size_t Iter = 0; Iter < IntervalIters; Iter++) {
    size_t ProgressSteps = 0;
    size_t SortedSteps = 0;

    std::fill(NeedUpdateNext.begin(), NeedUpdateNext.end(), false);

    uint32_t LBMLast = NOT_SET32;
    uint32_t PrevEnd = NOT_SET32;
    for (uint32_t Start = 0; Start + 1 < N; Start += StepSize) {
      const uint32_t End = std::min(N, Start + IntervalSize);
      if (End == PrevEnd)
        break;
      PrevEnd = End;

      if (!std::any_of(NeedUpdate.begin() + Start, NeedUpdate.begin() + End, [](bool x) { return x; }))
        continue;
      SortedSteps++;

      // Make sure LB values are pre-computed
      if (LBMLast == NOT_SET32 || End > LBMLast) {
        LBMLast = std::min(Start + Config.MaxDocsLB, N);
      }
      initLB(Documents.begin() + Start, Documents.begin() + LBMLast);

      uint64_t Delta = randomSearchIncremental(Documents.begin() + Start, Documents.begin() + End, RandomAttempts, true);
      if (Delta == 0)
        continue;

      SavedCrossings += Delta;
      ProgressSteps++;
      std::fill(NeedUpdateNext.begin() + Start, NeedUpdateNext.begin() + End, true);
    }
    NeedUpdate = NeedUpdateNext;

    LOG_IF(verbose >= 2, "    iteration %d/%d; progress steps: %d / %d", Iter + 1, IntervalIters, ProgressSteps, SortedSteps);
    if (ProgressSteps == 0)
      break;
  }

  if (Documents.size() <= Config.MaxDocsLB && SavedCrossings > 0) {
    SavedCrossings += tuneInterval(Documents.begin(), Documents.end(), true, OptLevelT::O1, 0);
  }

  CHECK(bucketsCorrect(Documents.begin(), Documents.end()));

  auto endTime = std::chrono::steady_clock::now();
  auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - beginTime).count();
  LOG_IF(verbose >= 1, "  removed %d crossings by tuning %d-intervals; running time is %d seconds", SavedCrossings, IntervalSize, durationSec);

  if (SavedCrossings > 0) {
    progressFunc("bp_" + std::to_string(seed) + "_tune" + std::to_string(IntervalSize) + "_int");
  }
  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::randomSearchWithSwaps(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t NumSwaps,
    const uint32_t MaxIters,
    const uint32_t MaxItersWithoutProgress) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  if (N > Config.MaxDocsLB)
    return 0;
  initLB(DocumentBegin, DocumentEnd);
  CHECK(0 < N && N <= Config.MaxDocsLB);

  auto beginTime = std::chrono::steady_clock::now();
  LOG_IF(verbose, "random search with %d swaps for %d documents", NumSwaps, N);

  // Apply swaps
  auto genOrderFunc = [&](
      const std::vector<Document*>::iterator& DocBegin,
      const std::vector<Document*>::iterator& DocEnd,
      const uint32_t Iter) {
    uint32_t MinBucket = DocBegin[0]->Bucket;
    for (auto It = DocBegin; It != DocEnd; It++)
      MinBucket = std::min(MinBucket, (*It)->Bucket);

    if (Iter == 0) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::AVERAGE);
      return;
    } else if (Iter == 1) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::MEDIAN);
      return;
    } else if (Iter == 2) {
      order(DocBegin, DocEnd, MinBucket, OrderAlgT::INPUT);
      return;
    }

    const uint32_t NumDocument = std::distance(DocBegin, DocEnd);
    for (uint32_t J = 0; J < NumSwaps; J++) {
      uint32_t D1 = Rand::next(NumDocument, RNG);
      uint32_t D2 = Rand::next(NumDocument, RNG);
      if (D1 == D2)
        continue;
      std::swap(DocBegin[D1]->Bucket, DocBegin[D2]->Bucket);
    }
    std::sort(DocBegin, DocEnd, Document::BucketOrder());
  };
  // Convergence
  auto convergeFunc = [&](const std::vector<uint64_t>& IterObjective) {
    return attemptsConverged(
      IterObjective,
      MaxIters,
      MaxItersWithoutProgress
    );
  };

  uint64_t SavedCrossings = randomSearch(DocumentBegin, DocumentEnd, true, genOrderFunc, convergeFunc);

  if (SavedCrossings > 0) {
    SavedCrossings += tuneInterval(DocumentBegin, DocumentEnd, true, OptLevelT::O1, 0);
  }

  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));

  auto endTime = std::chrono::steady_clock::now();
  auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - beginTime).count();
  LOG_IF(verbose >= 1, "  removed %d crossings by random %d-swaps; running time is %d seconds",
      SavedCrossings, NumSwaps, durationSec);

  if (SavedCrossings > 0) {
    progressFunc("bp_" + std::to_string(seed) + "_swaps" + std::to_string(NumSwaps));
  }
  return SavedCrossings;
}

template <class T>
std::pair<uint64_t, uint64_t> BalancedPartitioning<T>::computeCurOptCrossings(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  CHECK(NumDocuments <= Config.MaxDocsLB);
  // We assume the LBM values are initialized already!
  uint64_t CurCrossings = 0;
  uint64_t OptCrossings = 0;
  for (uint32_t I = 0; I < NumDocuments; I++) {
    CurCrossings += DocumentBegin[I]->getSelfCrossings();
    OptCrossings += DocumentBegin[I]->getSelfCrossings();
    for (uint32_t J = I + 1; J < NumDocuments; J++) {
      const uint64_t LB_LR = LBData[DocumentBegin[I]->LocalIndex][DocumentBegin[J]->LocalIndex];
      const uint64_t LB_RL = LBData[DocumentBegin[J]->LocalIndex][DocumentBegin[I]->LocalIndex];

      CurCrossings += LB_LR;
      OptCrossings += std::min(LB_LR, LB_RL);
    }
  }
  return {CurCrossings, OptCrossings};
}

template <class T>
bool BalancedPartitioning<T>::isOptimal(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  CHECK(NumDocuments <= Config.MaxDocsLB);
  // We assume the LBM values are initialized already!
  for (uint32_t I = 0; I < NumDocuments; I++) {
    for (uint32_t J = I + 1; J < NumDocuments; J++) {
      const uint64_t LB_LR = LBData[DocumentBegin[I]->LocalIndex][DocumentBegin[J]->LocalIndex];
      const uint64_t LB_RL = LBData[DocumentBegin[J]->LocalIndex][DocumentBegin[I]->LocalIndex];

      if (LB_LR > LB_RL)
        return false;
    }
  }
  return true;
}

template <class T>
uint64_t BalancedPartitioning<T>::sortIntervalOpt(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    bool IsIncremental) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N <= MaxDocsDP, "N = %d, MaxDocsDP = %d", N, MaxDocsDP);
  if (N <= 1)
    return 0;
  if (N == 2) {
    const T LB_01 = LBData[DocumentBegin[0]->LocalIndex][DocumentBegin[1]->LocalIndex];
    const T LB_10 = LBData[DocumentBegin[1]->LocalIndex][DocumentBegin[0]->LocalIndex];

    if (LB_01 <= LB_10)
      return 0;
    sortByOrder(DocumentBegin, DocumentEnd, {1, 0});
    return LB_01 - LB_10;
  }

  // Make sure LB values are pre-computed
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  // Stop early, if possible
  bool IsLastOptimal = true;
  bool IsFirstOptimal = true;
  for (uint32_t I = 1; I < N; I++) {
    if (IsFirstOptimal) {
      const T LB_0I = LBData[DocumentBegin[0]->LocalIndex][DocumentBegin[I]->LocalIndex];
      const T LB_I0 = LBData[DocumentBegin[I]->LocalIndex][DocumentBegin[0]->LocalIndex];
      if (LB_0I > LB_I0) {
        IsFirstOptimal = false;
      }
    }
    if (IsLastOptimal) {
      const T LB_IN = LBData[DocumentBegin[I - 1]->LocalIndex][DocumentBegin[N - 1]->LocalIndex];
      const T LB_NI = LBData[DocumentBegin[N - 1]->LocalIndex][DocumentBegin[I - 1]->LocalIndex];
      if (LB_IN > LB_NI) {
        IsLastOptimal = false;
      }
    }
  }
  // Stop early for incremental computations when the last vertex is optimal
  if (IsIncremental && IsLastOptimal)
    return 0;
  // No need to reorder the first one
  if (IsFirstOptimal) {
    return sortIntervalOpt(DocumentBegin + 1, DocumentEnd, IsIncremental);
  }
  // No need to reorder the last one
  if (IsLastOptimal) {
    CHECK(!IsIncremental);
    return sortIntervalOpt(DocumentBegin, DocumentEnd - 1, false);
  }

  uint32_t CurCrossings = 0;
  bool IsAllOptimal = true;

  for (uint32_t L = 0; L < N; L++) {
    LB_rev[L][L] = 0;
    for (uint32_t R = L + 1; R < N; R++) {
      const T LB_LR = LBData[DocumentBegin[L]->LocalIndex][DocumentBegin[R]->LocalIndex];
      const T LB_RL = LBData[DocumentBegin[R]->LocalIndex][DocumentBegin[L]->LocalIndex];
      const T LB_min = std::min(LB_LR, LB_RL);

      LB_rev[R][L] = LB_LR - LB_min;
      LB_rev[L][R] = LB_RL - LB_min;
      CHECK(LB_rev[R][L] == 0 || LB_rev[L][R] == 0);

      CurCrossings += (LB_LR - LB_min);
      if (LB_LR != LB_RL) {
        IsAllOptimal = false;
      }
    }
  }
  // Don't run it if we already have the optimum
  if (IsAllOptimal)
    return 0;

  // Trying to find an optimal spot for the last element
  if (IsIncremental) {
    for (uint32_t I = 0; I < N; I++) {
      const uint32_t L = N - 1 - I;
      // placing before L-th: [0 ... Last L ... N - 2]
      bool BestForLast = true;
      for (uint32_t J = 0; J < L; J++) {
        if (LB_rev[N - 1][J] > 0) {
          BestForLast = false;
          break;
        }
      }
      for (uint32_t J = L; J < N - 1; J++) {
        if (LB_rev[J][N - 1] > 0) {
          BestForLast = false;
          break;
        }
      }
      if (!BestForLast)
        continue;
      // Found a best position for the last element
      std::vector<uint32_t> Order;
      Order.reserve(N);
      for (uint32_t J = 0; J < L; J++) {
        Order.push_back(J);
        CHECK(LB_rev[N - 1][J] == 0);
      }
      Order.push_back(N - 1);
      uint32_t SavedCrossings = 0;
      for (uint32_t J = L; J < N - 1; J++) {
        Order.push_back(J);
        CHECK(LB_rev[J][N - 1] == 0);
        SavedCrossings += LB_rev[N - 1][J];
      }
      sortByOrder(DocumentBegin, DocumentEnd, Order);
      return SavedCrossings;
    }
  }

  using MaskType = uint32_t;
  constexpr MaskType ONE = static_cast<MaskType>(1);
  const MaskType MAX_MASK = ONE << N;
  const MaskType FIN_MASK = MAX_MASK - 1;

  // Re-compute DPPrevs if needed
  if (DPPrevs.size() < MAX_MASK) {
    LOG_IF(verbose >= 2, "Init DPPrevs for MAX_MASK = %d", MAX_MASK);
    DPPrevs.clear();
    DPPrevs.resize(MAX_MASK);
    for (MaskType Set = 1; Set < MAX_MASK; Set++) {
      auto& Prevs = DPPrevs[Set];
      Prevs.reserve(__builtin_popcount(Set));
      MaskType Mask = Set;
      while (Mask > 0) {
        Prevs.push_back(__builtin_ffs(Mask) - 1);
        Mask &= (Mask - 1);
      }
    }
  }
  // Pre-compute masks
  while (DPMasks.size() < N) {
    DPMasks.push_back({});
  }
  if (DPMasks[N - 1].empty()) {
    LOG_IF(verbose >= 2, "Init DPMasks for N - 1 = %d", N - 1);
    DPMasks[N - 1].reserve(10000);
    for (MaskType Set = 1; Set < MAX_MASK; Set++) {
      const auto& Prevs = DPPrevs[Set];
      if (Prevs.size() <= (N + 1) / 2) {
        DPMasks[N - 1].push_back(Set);
      }
    }
    DPMasks[N - 1].shrink_to_fit();
  }

  // Pre-compute successors and predecessors
  MaskType Pred[MaxDocsDP];
  std::fill(std::begin(Pred), std::begin(Pred) + MaxDocsDP, 0);
  MaskType Succ[MaxDocsDP];
  std::fill(std::begin(Succ), std::begin(Succ) + MaxDocsDP, 0);

  for (uint32_t L = 0; L < N; L++) {
    for (uint32_t R = 0; R < N; R++) {
      if (L == R) continue;
      if (LB_rev[R][L] <= LB_rev[L][R]) {
      //if (LB[L][R] <= LB[R][L]) {
        Pred[R] |= (static_cast<MaskType>(1) << L);
        Succ[L] |= (static_cast<MaskType>(1) << R);
      }
    }
  }

  // LastSum[I] == sum when I-th is the last
  uint32_t LastSum[MaxDocsDP];
  for (uint32_t I = 0; I < N; I++) {
    LastSum[I] = 0;
    for (uint32_t J = 0; J < N; J++) {
      LastSum[I] += LB_rev[I][J];
    }
  }

  // Initialize the DP entries
  std::fill(DP.begin(), DP.begin() + MAX_MASK, MAX_32);
  DP[0] = 0;

  // Walk through solutions using a bitmask to represent the state
  const auto& Masks = DPMasks[N - 1];
  const uint32_t NumMasks = Masks.size();
  for (uint32_t SetIdx = 0; SetIdx < NumMasks; SetIdx++) {
    const MaskType Set = Masks[SetIdx];
    // Extract set elements.
    const auto& Prevs = DPPrevs[Set];
    const uint16_t NumPrevs = Prevs.size();
    // if (NumPrevs > (N + 1) / 2)
    //   continue;
    CHECK(NumPrevs <= (N + 1) / 2);

    // Try to find a forced Last
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t Last = Prevs[NumPrevs - LastIdx - 1];
      const MaskType OtherSet = Set ^ (ONE << Last);
      if ((OtherSet & Pred[Last]) == OtherSet) {
        DPFirst[Set] = NOT_SET16;
        DPLast[Set] = Last;
        DP[Set] = DP[OtherSet];
        break;
      }
    }
    if (DP[Set] != MAX_32)
      continue;

    // Try to find a forced First
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t First = Prevs[LastIdx];
      const MaskType OtherSet = Set ^ (ONE << First);
      if ((OtherSet & Succ[First]) == OtherSet) {
        DPFirst[Set] = First;
        DPLast[Set] = NOT_SET16;
        DP[Set] = DP[OtherSet];
        break;
      }
    }
    if (DP[Set] != MAX_32)
      continue;

    // Traverse each possibility of Last document visited in this layout
    uint32_t BestSetCrossings = MAX_32;
    uint16_t BestSetLast = 0;
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t Last = Prevs[NumPrevs - LastIdx - 1];
      const MaskType PrevSet = Set ^ (ONE << Last);

      uint32_t LastCrossings = 0;

      uint16_t PrevIdx = 0;
      while (PrevIdx < NumPrevs) {
        const uint16_t Prev = Prevs[PrevIdx++];
        LastCrossings += LB_rev[Last][Prev];
        //__builtin_prefetch(&LB_rev[Last][Prev + 1], 0, 3);
      }

      const uint32_t NumCrossings = DP[PrevSet] + LastCrossings;
      if (BestSetCrossings > NumCrossings) {
        BestSetCrossings = NumCrossings;
        BestSetLast = Last;
        if (BestSetCrossings == 1)
          break;
      }
    }

    DP[Set] = BestSetCrossings;
    DPLast[Set] = BestSetLast;
    DPFirst[Set] = NOT_SET16;

    // Couldn't find a better solution; stop early
    const MaskType Set2 = FIN_MASK ^ Set;
    if ((DP[Set] >= CurCrossings) || (DP[Set2] != MAX_32 && DP[Set] + DP[Set2] >= CurCrossings)) {
      return 0;
    }
  }

  // Calculate the best value
  uint32_t BestCrossings = MAX_32;
  const uint32_t K = (N + 1) / 2;
  const MaskType InitSet = (ONE << K) - 1;
  const MaskType LastSet = (ONE << N) - (ONE << (N - K));
  for (MaskType Set = InitSet, r, c; Set <= LastSet; c = Set & -Set, r = Set + c, Set = r | (((r ^ Set) >> 2) / c)) {
    const MaskType Set2 = FIN_MASK ^ Set;
    CHECK(DP[Set] < MAX_32 && DP[Set2] < MAX_32);

    uint32_t AllCrossings = DP[Set] + DP[Set2];
    if (AllCrossings >= CurCrossings)
      continue;

    const auto& Prevs = DPPrevs[Set];
    const uint16_t NumPrevs = Prevs.size();
    const auto& Nexts = DPPrevs[Set2];
    const uint16_t NumNexts = Nexts.size();

    for (uint16_t NextIdx = 0; NextIdx < NumNexts; NextIdx++) {
      const uint16_t Next = Nexts[NextIdx];
      for (uint16_t PrevIdx = 0; PrevIdx < NumPrevs; PrevIdx++) {
        const uint16_t Prev = Prevs[PrevIdx];
        AllCrossings += LB_rev[Next][Prev];
      }
    }

    if (BestCrossings > AllCrossings) {
      BestCrossings = AllCrossings;
    }
  }

  // Stop early, if no improvement
  if (BestCrossings >= CurCrossings) {
    return 0;
  }

  // Continue with the remaining computation
  for (MaskType Set = InitSet; Set < MAX_MASK; Set++) {
    // Extract set elements.
    const auto& Prevs = DPPrevs[Set];
    const uint16_t NumPrevs = Prevs.size();
    if (NumPrevs <= (N + 1) / 2)
      continue;

    // Try to find a forced Last
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t Last = Prevs[NumPrevs - LastIdx - 1];
      const MaskType OtherSet = Set ^ (static_cast<MaskType>(1) << Last);
      if ((OtherSet & Pred[Last]) == OtherSet) {
        DPFirst[Set] = NOT_SET16;
        DPLast[Set] = Last;
        DP[Set] = DP[OtherSet];
        break;
      }
    }
    if (DP[Set] != MAX_32)
      continue;

    // Try to find a forced First
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t First = Prevs[LastIdx];
      const MaskType OtherSet = Set ^ (static_cast<MaskType>(1) << First);
      if ((OtherSet & Succ[First]) == OtherSet) {
        DPFirst[Set] = First;
        DPLast[Set] = NOT_SET16;
        DP[Set] = DP[OtherSet];
        break;
      }
    }
    if (DP[Set] != MAX_32)
      continue;

    const MaskType Set2 = FIN_MASK ^ Set;
    const auto& Nexts = DPPrevs[Set2];
    const uint16_t NumNexts = Nexts.size();

    // Traverse each possibility of Last document visited in this layout
    uint32_t BestSetCrossings = MAX_32;
    uint16_t BestSetLast = 0;
    for (uint16_t LastIdx = 0; LastIdx < NumPrevs; LastIdx++) {
      const uint16_t Last = Prevs[NumPrevs - LastIdx - 1];
      const MaskType PrevSet = Set ^ (ONE << Last);

      // Iterating over remaining items, Next, and compute LastSum[Last] - sum(LB[Next][Last])
      uint32_t LastCrossings = 0;

      uint16_t NextIdx = 0;
      while (NextIdx < NumNexts) {
        const uint16_t Next = Nexts[NextIdx++];
        LastCrossings += LB_rev[Last][Next];
      }
      // CHECK(LastSum[Last] >= LastCrossings, "LastSum[Last] = %d; LastCrossings = %d", LastSum[Last], LastCrossings);
      LastCrossings = LastSum[Last] - LastCrossings;

      const uint32_t NumCrossings = DP[PrevSet] + LastCrossings;
      if (BestSetCrossings > NumCrossings) {
        BestSetCrossings = NumCrossings;
        BestSetLast = Last;
        if (BestSetCrossings == 1)
          break;
      }
    }

    DP[Set] = BestSetCrossings;
    DPLast[Set] = BestSetLast;
    DPFirst[Set] = NOT_SET16;
  }
  CHECK(CurCrossings > DP[FIN_MASK] && DP[FIN_MASK] == BestCrossings);

  CHECK(CurCrossings >= DP[FIN_MASK]);
  const uint64_t SavedCrossings = CurCrossings - DP[FIN_MASK];

  if (SavedCrossings == 0)
    return 0;

  // Extract the permutation
  std::vector<uint32_t> BestPermL;
  BestPermL.reserve(N);
  std::vector<uint32_t> BestPermF;
  BestPermF.reserve(N);
  MaskType CurSet = FIN_MASK;
  while (CurSet > 0) {
    CHECK(DPLast[CurSet] == NOT_SET16 || DPFirst[CurSet] == NOT_SET16);
    if (DPLast[CurSet] != NOT_SET16) {
      BestPermL.push_back(DPLast[CurSet]);
      CurSet ^= (1 << BestPermL.back());
    } else {
      CHECK(DPFirst[CurSet] != NOT_SET16);
      BestPermF.push_back(DPFirst[CurSet]);
      CurSet ^= (1 << BestPermF.back());
    }
  }
  std::reverse(BestPermL.begin(), BestPermL.end());

  std::vector<uint32_t> BestPerm;
  BestPerm.insert(BestPerm.end(), BestPermF.begin(), BestPermF.end());
  BestPerm.insert(BestPerm.end(), BestPermL.begin(), BestPermL.end());
  CHECK(BestPerm.size() == N);

  // Sort by the permutation
  sortByOrder(DocumentBegin, DocumentEnd, BestPerm);
  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::sortIntervalsOpt(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<void(const std::string& name)> progressFunc,
    const uint32_t IntervalSize,
    const uint32_t IntervalIter,
    const OptLevelT OptLevel,
    const size_t verbose) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N > 0);
  CHECK(IntervalSize <= Config.MaxDocsLB);
  LOG_IF(verbose, "sorting intervals with size = %d for %d docs", IntervalSize, N);

  uint64_t SavedCrossings = 0;

  std::vector<bool> NeedUpdate(N, true);
  std::vector<bool> NeedUpdateNext(N, false);
  size_t StepsFromProgress = 0;
  for (size_t Iter = 0; Iter < IntervalIter; Iter++) {
    size_t ProgressSteps = 0;
    size_t TotalSteps = 0;
    size_t SortedSteps = 0;
    std::fill(NeedUpdateNext.begin(), NeedUpdateNext.end(), false);
    uint32_t LBMLast = NOT_SET32;
    bool IsIncremental = false;
    uint32_t PrevEnd = NOT_SET32;
    for (uint32_t Start = 0; Start + 1 < N; Start++) {
      const uint32_t End = std::min(N, Start + IntervalSize);
      if (End == PrevEnd)
        break;
      PrevEnd = End;

      TotalSteps++;
      if (!std::any_of(NeedUpdate.begin() + Start, NeedUpdate.begin() + End, [](bool x) { return x; }))
        continue;
      SortedSteps++;

      // Make sure LB values are pre-computed
      if (LBMLast == NOT_SET32 || End > LBMLast) {
        LBMLast = std::min(Start + Config.MaxDocsLB, N);
      }
      initLB(DocumentBegin + Start, DocumentBegin + LBMLast);

      uint64_t Delta = sortIntervalOpt(DocumentBegin + Start, DocumentBegin + End, IsIncremental);
      IsIncremental = true;
      if (Delta == 0) {
        continue;
      }

      SavedCrossings += Delta;
      ProgressSteps++;
      std::fill(NeedUpdateNext.begin() + Start, NeedUpdateNext.begin() + End, true);
    }
    NeedUpdate = NeedUpdateNext;

    LOG_IF(verbose, "    iteration %d/%d; progress steps: %d / %d", Iter + 1, IntervalIter, ProgressSteps, SortedSteps);
    if (ProgressSteps == 0)
      break;
    StepsFromProgress += ProgressSteps;
    if (Iter + 1 < IntervalIter && StepsFromProgress * 8 >= TotalSteps) {
      progressFunc("bp_" + std::to_string(seed) + "_pp_" + std::to_string(Config.PostIntervalSize));
      StepsFromProgress = 0;
    }
  }
  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::mergeTwoIntervals(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    std::function<bool(uint32_t, const Document*)> IsFirstFilter) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N <= Config.MaxDocsLB);
  LOG_IF(verbose >= 3, "mergeTwoIntervals for %d docs", N);
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  std::vector<uint32_t> Order1;
  Order1.reserve(N);
  std::vector<uint32_t> Order2;
  Order2.reserve(N);
  for (uint32_t DocIdx = 0; DocIdx < N; DocIdx++) {
    if (IsFirstFilter(DocIdx, DocumentBegin[DocIdx])) {
      Order1.push_back(DocIdx);
    } else {
      Order2.push_back(DocIdx);
    }
  }
  const uint32_t N1 = Order1.size();
  const uint32_t N2 = Order2.size();
  if (N1 == 0 || N2 == 0)
    return 0;
  CHECK(N1 > 0 && N2 > 0 && N == N1 + N2);

  // Allocate the data and compute partial sums
  allocateMatrix(LBSum1, N1 + 1, N2 + 1, Config.MaxDocsLB + 1, 0);
  allocateMatrix(LBSum2, N1 + 1, N2 + 1, Config.MaxDocsLB + 1, 0);
  uint64_t CurCrossings = 0;
  uint64_t OptCrossings = 0;
  for (uint32_t I = 0; I < N1; I++) {
    for (uint32_t J = 0; J < N2; J++) {
      const uint32_t L = Order1[I];
      const uint32_t R = Order2[J];

      const uint64_t LB_LR = LBData[DocumentBegin[L]->LocalIndex][DocumentBegin[R]->LocalIndex];
      const uint64_t LB_RL = LBData[DocumentBegin[R]->LocalIndex][DocumentBegin[L]->LocalIndex];

      LBSum1[I + 1][J + 1] = LBSum1[I + 1][J] + LB_RL;
      LBSum2[I + 1][J + 1] = LBSum2[I][J + 1] + LB_LR;

      CurCrossings += L < R ? LB_LR : LB_RL;
      OptCrossings += std::min(LB_LR, LB_RL);
    }
  }

  // Don't run it if we already have the optimum
  if (CurCrossings == OptCrossings)
    return 0;

  // DP[i][j] == min crossings after merging (Order1[0] .. Order1[i-1]) with (Order2[0] .. Order2[j-1])
  allocateMatrix(MergeDP, N1 + 1, N2 + 1, Config.MaxDocsLB + 1, NOT_SET64);

  // Init
  MergeDP[0][0] = 0;
  for (uint32_t I = 1; I <= N1; I++) {
    MergeDP[I][0] = 0;
  }
  for (uint32_t J = 1; J <= N2; J++) {
    MergeDP[0][J] = 0;
  }

  // Apply DP
  for (uint32_t I = 1; I <= N1; I++) {
    for (uint32_t J = 1; J <= N2; J++) {
      CHECK(MergeDP[I][J] == NOT_SET64 && MergeDP[I - 1][J] != NOT_SET64 && MergeDP[I][J - 1] != NOT_SET64);

      // Place Order1[I-1] last in the order
      const uint64_t CrossLastI = MergeDP[I - 1][J] + LBSum1[I][J];
      // Place Order2[J] last in the order
      const uint64_t CrossLastJ = MergeDP[I][J - 1] + LBSum2[I][J];

      MergeDP[I][J] = std::min(CrossLastI, CrossLastJ);
    }
  }

  const uint64_t BestCrossing = MergeDP[N1][N2];
  CHECK(BestCrossing <= CurCrossings);
  // Stop early if no better solution is found
  if (BestCrossing == CurCrossings)
    return 0;

  // Extract the permutation
  uint32_t NowI = N1;
  uint32_t NowJ = N2;
  std::vector<uint32_t> BestOrder;
  BestOrder.reserve(N);
  while (NowI > 0 || NowJ > 0) {
    if (NowI > 0 && MergeDP[NowI][NowJ] == MergeDP[NowI - 1][NowJ] + LBSum1[NowI][NowJ]) {
      NowI--;
      BestOrder.push_back(Order1[NowI]);
    } else {
      CHECK(NowJ > 0 && MergeDP[NowI][NowJ] == MergeDP[NowI][NowJ - 1] + LBSum2[NowI][NowJ]);
      NowJ--;
      BestOrder.push_back(Order2[NowJ]);
    }
  }
  std::reverse(BestOrder.begin(), BestOrder.end());
  CHECK(BestOrder.size() == N);

  // Sort with respect to the best order
  sortByOrder(DocumentBegin, DocumentEnd, BestOrder);

  LOG_IF(verbose >= 2, "  eliminated %d (%d -> %d) crossings by merging two intervals", CurCrossings - BestCrossing, CurCrossings, BestCrossing);
  return CurCrossings - BestCrossing;
}

template <class T>
uint64_t BalancedPartitioning<T>::mergeTwoIntervals(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentMid,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  const uint32_t N1 = std::distance(DocumentBegin, DocumentMid);
  const uint32_t N2 = std::distance(DocumentMid, DocumentEnd);
  if (N1 == 0 || N2 == 0)
    return 0;

  std::function<bool(uint32_t, const Document*)> N1Filter = [&](uint32_t DocIdx, const Document* Doc) {
    CHECK(DocIdx < N1 + N2);
    return DocIdx < N1;
  };
  return mergeTwoIntervals(DocumentBegin, DocumentEnd, N1Filter);
}

template <class T>
uint64_t BalancedPartitioning<T>::mergeTwoIntervalsByDegree(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    uint32_t MaxDegree) const {
  std::function<bool(uint32_t, const Document*)> Filter = [&](uint32_t DocIdx, const Document* Doc) {
    return Doc->degree() <= MaxDegree;
  };
  return mergeTwoIntervals(DocumentBegin, DocumentEnd, Filter);
}

template <class T>
uint64_t BalancedPartitioning<T>::mergeTwoIntervalsOddEven(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  std::function<bool(uint32_t, const Document*)> Filter = [&](uint32_t DocIdx, const Document* Doc) {
    return DocIdx % 2 == 0;
  };
  return mergeTwoIntervals(DocumentBegin, DocumentEnd, Filter);
}

template <class T>
uint64_t BalancedPartitioning<T>::mergeIntervalPairs(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t MaxDiffSize) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  if (N > Config.MaxDocsLB)
    return 0;

  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  // Don't run it if we already have the optimum
  if (isOptimal(DocumentBegin, DocumentEnd))
    return 0;

  uint64_t SavedCrossings = 0;
  const uint32_t StepSize = 1;
  for (uint32_t I = 1; I < N; I += StepSize) {
    const uint32_t BI = I >= MaxDiffSize ? I - MaxDiffSize : 0;
    const uint32_t EI = std::min(I + MaxDiffSize, N);
    CHECK(BI < I && I < EI);

    const std::vector<Document*>::iterator DocBegin = DocumentBegin + BI;
    const std::vector<Document*>::iterator DocMid = DocumentBegin + I;
    const std::vector<Document*>::iterator DocEnd = DocumentBegin + EI;
    SavedCrossings += mergeTwoIntervals(DocBegin, DocMid, DocEnd);
  }

  return SavedCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::adjustSingle(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t FirstIndex,
    std::vector<uint64_t>& SumL,
    std::vector<uint64_t>& SumR) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N > 1);

  const uint32_t FirstLI = DocumentBegin[FirstIndex]->LocalIndex;

  // Check if already optimal
  bool IsOptimal = true;
  for (uint32_t I = 0; I < FirstIndex; I++) {
    const uint32_t OtherLI = DocumentBegin[I]->LocalIndex;
    if (LBData[OtherLI][FirstLI] > LBData[FirstLI][OtherLI]) {
      IsOptimal = false;
      break;
    }
  }
  if (IsOptimal) {
    for (uint32_t I = FirstIndex + 1; I < N; I++) {
      const uint32_t OtherLI = DocumentBegin[I]->LocalIndex;
      if (LBData[FirstLI][OtherLI] > LBData[OtherLI][FirstLI]) {
        IsOptimal = false;
        break;
      }
    }
  }
  if (IsOptimal) {
    return 0;
  }

  // Apply the optimization
  const uint32_t N2 = N - 1;

  // SumL[I] = sum LB(J, Doc) for J <= I
  SumL[0] = 0;
  for (uint32_t I = 0; I < FirstIndex; I++) {
    // I is the insertion index
    SumL[I + 1] = SumL[I] + LBData[DocumentBegin[I]->LocalIndex][FirstLI];
  }
  for (uint32_t I = FirstIndex; I < N2; I++) {
    SumL[I + 1] = SumL[I] + LBData[DocumentBegin[I + 1]->LocalIndex][FirstLI];
  }

  // SumR[I] = sum LB(J, Doc) for J >= I
  SumR[N2] = 0;
  for (uint32_t I = 0; I + FirstIndex < N2; I++) {
    const uint32_t J = N2 - 1 - I;
    SumR[J] = SumR[J + 1] + LBData[FirstLI][DocumentBegin[J + 1]->LocalIndex];
  }
  for (uint32_t I = N2 - FirstIndex; I < N2; I++) {
    const uint32_t J = N2 - 1 - I;
    SumR[J] = SumR[J + 1] + LBData[FirstLI][DocumentBegin[J]->LocalIndex];
  }

  const uint64_t CurCrossings = SumL[FirstIndex] + SumR[FirstIndex];
  uint64_t BestCrossings = MAX_64;
  uint32_t BestPos = NOT_SET32;
  for (uint32_t I = 0; I < N; I++) {
    // insertion point is right before I-th in Order2
    const uint64_t CrossingsI = SumL[I] + SumR[I];
    if (BestCrossings > CrossingsI) {
      BestCrossings = CrossingsI;
      BestPos = I;
    }
  }
  CHECK(BestCrossings <= CurCrossings);

  // Stop early if no better solution is found
  if (BestCrossings == CurCrossings)
    return 0;
  CHECK(BestPos != FirstIndex);

  // Insert the element at position BestPos
  const uint32_t FirstBucket = DocumentBegin[FirstIndex]->Bucket;
  const uint32_t BestBucket = DocumentBegin[BestPos]->Bucket;
  if (FirstIndex + 1 == BestPos) {
    std::swap(DocumentBegin[FirstIndex], DocumentBegin[FirstIndex + 1]);
    DocumentBegin[FirstIndex]->Bucket = FirstBucket;
    DocumentBegin[FirstIndex + 1]->Bucket = FirstBucket + 1;
  } else if (BestPos + 1 == FirstIndex) {
    CHECK(FirstBucket == BestBucket + 1);
    std::swap(DocumentBegin[BestPos], DocumentBegin[BestPos + 1]);
    DocumentBegin[BestPos]->Bucket = BestBucket;
    DocumentBegin[BestPos + 1]->Bucket = BestBucket + 1;
  } else if (FirstIndex < BestPos) {
    std::rotate(DocumentBegin + FirstIndex, DocumentBegin + FirstIndex + 1, DocumentBegin + BestPos + 1);
    for (uint32_t I = FirstIndex, J = 0; I <= BestPos; I++, J++) {
      DocumentBegin[I]->Bucket = FirstBucket + J;
    }
  } else {
    CHECK(FirstIndex > BestPos);
    std::rotate(DocumentBegin + BestPos, DocumentBegin + FirstIndex, DocumentBegin + FirstIndex + 1);
    for (uint32_t I = BestPos, J = 0; I <= FirstIndex; I++, J++) {
      DocumentBegin[I]->Bucket = BestBucket + J;
    }
  }

  LOG_IF(verbose >= 4, "  eliminated %d (%d -> %d) crossings by reordering single document",
      CurCrossings - BestCrossings, CurCrossings, BestCrossings);
  return CurCrossings - BestCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::adjustSingle(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const OptLevelT OptLevel) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  CHECK(0 < NumDocuments && NumDocuments <= Config.MaxDocsLB);
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));
  LOG_IF(verbose >= 3, "adjustSingle for %d docs", NumDocuments);

  const auto [CurCrossings, OptCrossings] = computeCurOptCrossings(DocumentBegin, DocumentEnd);
  // Don't run it if we already have the optimum
  if (CurCrossings == OptCrossings)
    return 0;
  uint64_t NewCrossings = CurCrossings;

  // Pre-allocate data for the algorithm
  std::vector<uint64_t> SumL(NumDocuments);
  std::vector<uint64_t> SumR(NumDocuments);

  size_t NumUpdateIters = 0;
  while (true) {
    size_t NumUpdates = 0;
    for (uint32_t I = 0; I < NumDocuments; I++) {
      uint64_t Delta = adjustSingle(DocumentBegin, DocumentEnd, I, SumL, SumR);
      if (Delta > 0) {
        NewCrossings -= Delta;
        NumUpdates++;
        I--;
      }
    }
    if (NumUpdates == 0 || NewCrossings == OptCrossings)
      break;
    NumUpdateIters++;

    if (OptLevel <= OptLevelT::O1 && Config.PostTuneInt == 0 && NumUpdateIters >= 3) {
      break;
    }
    if (OptLevel <= OptLevelT::O1 && Config.PostTuneInt == 1 && NumUpdateIters >= 9) {
      break;
    }
  }

  CHECK(NewCrossings <= CurCrossings);
  CHECK(NumUpdateIters == 0 || NewCrossings < CurCrossings);

  if (NewCrossings < CurCrossings) {
    LOG_IF(verbose >= 2, "  eliminated %d (%'lld -> %'lld) crossings by reordering single documents in %d iterations",
        CurCrossings - NewCrossings, CurCrossings, NewCrossings, NumUpdateIters);
  }
  return CurCrossings - NewCrossings;
}

template class BalancedPartitioning<uint8_t>;
template class BalancedPartitioning<uint16_t>;
template class BalancedPartitioning<uint32_t>;
template class BalancedPartitioning<uint64_t>;
