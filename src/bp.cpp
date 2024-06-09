#include "bp.h"
#include "logging.h"
#include "common.h"

void BalancedPartitioningConfig::setOrderAlg(const std::string& alg) {
  if (alg == "input")
    OrderAlg = OrderAlgT::INPUT;
  else if (alg == "median")
    OrderAlg = OrderAlgT::MEDIAN;
  else if (alg == "average")
    OrderAlg = OrderAlgT::AVERAGE;
  else if (alg == "local-median")
    OrderAlg = OrderAlgT::LOCAL_MEDIAN;
  else if (alg == "local-average")
    OrderAlg = OrderAlgT::LOCAL_AVERAGE;
  else
    ERROR("unknown order-alg: " + alg);
}

void BalancedPartitioningConfig::setSplitAlg(const std::string& alg) {
  if (alg == "input")
    SplitAlg = SplitAlgT::INPUT;
  else if (alg == "random")
    SplitAlg = SplitAlgT::RANDOM;
  else if (alg == "median")
    SplitAlg = SplitAlgT::MEDIAN;
  else if (alg == "average")
    SplitAlg = SplitAlgT::AVERAGE;
  else
    ERROR("unknown split-alg: " + alg);
}

template <class T>
BalancedPartitioning<T>::BalancedPartitioning(const BalancedPartitioningConfig& Config,
                                              const bool UseDenseLBData)
    : Config(Config), UseDenseLBData(UseDenseLBData) {}

template <class T>
void BalancedPartitioning<T>::initialize(std::vector<Document *>& Documents) {
  MaxUtility = NOT_SET32;
  for (const Document* Doc : Documents) {
    CHECK(Doc->degree() > 0);
    const uint32_t MaxVec = Doc->getOrigEdges().back().utility;
    if (MaxUtility == NOT_SET32 || MaxUtility < MaxVec)
      MaxUtility = MaxVec;
  }
  CHECK(MaxUtility != NOT_SET32);

  stree.resize(2 * (MaxUtility + 1));
  DP.resize(MaxDPSize);
  DPLast.resize(MaxDPSize);
  DPFirst.resize(MaxDPSize);

  RNG.seed(seed);

  std::vector<uint32_t> UNodes;
  for (uint32_t I = 0; I < Documents.size(); I++) {
    Document* Doc = Documents[I];
    Doc->Bucket = NOT_SET32;
    UNodes.clear();
    UNodes.reserve(Doc->degree());
    for (const auto [U, _] : Doc->getEdges()) {
      UNodes.push_back(U);
    }
    Doc->AllMedian = median(UNodes);
    Doc->AllAverage = average(UNodes);
  }

  // Init lower bounds
  if (Documents.size() <= Config.MaxDocsLB && LowerBound == NOT_SET64) {
    const auto [_, OptCrossings] = computeCurOptCrossings(Documents.begin(), Documents.end());
    LowerBound = OptCrossings;
    LOG_IF(verbose >= 1, " original lower bound: %'lld", LowerBound);
  }
}

template <class T>
void BalancedPartitioning<T>::run(std::vector<Document *>& Documents) {
  auto beginTime = std::chrono::steady_clock::now();

  initialize(Documents);

  LOG_IF(verbose, "Partitioning %d documents with %d utilities using depth %d, %d iterations per split, and %d split attempts",
    Documents.size(), MaxUtility + 1, Config.SplitDepth, Config.IterationsPerSplit, Config.SplitAttempts);

  bisect(begin(Documents), end(Documents), 0, 1, 0, true);
  // Make sure all utilities are back
  restoreOrigUtilities(begin(Documents), end(Documents));
  // Make sure the final buckets are correct
  for (uint32_t I = 0; I < Documents.size(); I++) {
    CHECK(Documents[I]->Bucket == I);
  }

  auto endTime = std::chrono::steady_clock::now();
  auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - beginTime).count();
  LOG_IF(verbose, "Balanced partitioning completed; running time is %d seconds", durationSec);
}

template <class T>
void BalancedPartitioning<T>::bisect(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    const uint32_t RecDepth, const uint32_t RootBucket,
    const uint32_t Offset, const bool SplitBySCCs) {
  LOG_IF(verbose >= 2, "BalancedPartitioning<T>::bisect with RecDepth = %d; RootBucket = %d; Offset = %d; SplitBySCCs = %d; N = %d",
    RecDepth, RootBucket, Offset, SplitBySCCs, std::distance(DocumentBegin, DocumentEnd));

  // No work to do
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  if (NumDocuments == 0)
    return;

  // Reached the lowest level of the recursion tree
  if (RecDepth >= Config.SplitDepth || NumDocuments <= 1 || NumDocuments <= Config.LeafInterval) {
    order(DocumentBegin, DocumentEnd, Offset, Config.OrderAlg);
    return;
  }

  const bool UseStrongMoveGains = NumDocuments <= Config.MaxDocsLB;

  // Initialize LB matrix to use throughout the bisection step
  if (UseStrongMoveGains) {
    initLB(DocumentBegin, DocumentEnd);
  }

  // Stop early if we already have the optimum
  if (UseStrongMoveGains && isOptimal(DocumentBegin, DocumentEnd)) {
    order(DocumentBegin, DocumentEnd, Offset, OrderAlgT::NONE);
    return;
  }

  // Split the graph into SCCs and apply reordering separately
  if (UseStrongMoveGains && SplitBySCCs) {
    std::vector<IterRange> Ranges;
    splitAndReorderBySCCs(DocumentBegin, DocumentEnd, Ranges);
    CHECK(Ranges.size() >= 1);
    if (Ranges.size() > 1) {
      uint32_t NumProcessed = 0;
      for (const auto [DocBegin, DocEnd] : Ranges) {
        bisect(DocBegin, DocEnd, RecDepth, RootBucket, Offset + NumProcessed, false);
        NumProcessed += std::distance(DocBegin, DocEnd);
      }
      CHECK(NumProcessed == NumDocuments);
      return;
    }
  }

  // Set variables
  this->UseStrongMoveGains = UseStrongMoveGains;
  const uint32_t LeftBucket = 2 * RootBucket;
  this->LeftBucket = LeftBucket;
  const uint32_t RightBucket = 2 * RootBucket + 1;
  this->RightBucket = RightBucket;

  //////////////////////////////////////////////////////////////////////////////
  // Start of Optimization
  //////////////////////////////////////////////////////////////////////////////

  // Adjust document adjacencies: renumber utilities and drop obsolete ones
  const uint32_t NumUtilities = updateUtilities(DocumentBegin, DocumentEnd);

  // The order of ids
  std::vector<uint32_t> OrgIdOrder;
  // Save original number of crossings
  auto OrgResult = orderForBacktrack(DocumentBegin, DocumentEnd, Offset, NumUtilities, OrgIdOrder);
  const uint64_t OrgNumCrossings = OrgResult.first;
  // Stop if there are no crossings
  if (OrgNumCrossings == 0) {
    order(DocumentBegin, DocumentEnd, Offset, OrgResult.second);
    return;
  }

  // Restore the order after backtracking
  std::sort(DocumentBegin, DocumentEnd, Document::LocalMedianOrder());
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    (*It)->Bucket = RootBucket;
  }

  LOG_IF(verbose >= 2, "  bisecting at depth %2d, bucket %3d: |D| = %4d, |U| = %4d, orig_crossings = %'lld [alg-%d]",
      RecDepth, RootBucket, NumDocuments, NumUtilities, OrgNumCrossings, OrgResult.second);


  // Define signatures
  SignaturesType Signatures(NumUtilities);

  uint64_t BestObjective = NOT_SET64;
  std::vector<uint64_t> IterObjective;
  // if the last several attempts produce the same result, stop early
  const uint32_t MaxEqualAttempts = 3;
  for (size_t Iter = 0; Iter < Config.SplitAttempts; Iter++) {
    // Initialize 2 buckets
    split(DocumentBegin, DocumentEnd, LeftBucket);

    // Initialize signatures
    initializeSignatures(Signatures, DocumentBegin, DocumentEnd);

    if (UseStrongMoveGains && Iter == 0) {
      initializeLBMatrix(DocumentBegin, DocumentEnd);
    }

    // Apply iterations to improve the objective
    runIterations(Signatures, DocumentBegin, DocumentEnd, RecDepth);

    if (Config.SplitAttempts == 1)
      break;

    const uint64_t CurObjective = objective(Signatures, DocumentBegin, DocumentEnd);
    if (CurObjective == NOT_SET64)
      break;

    LOG_IF(verbose >= 2, "    attempt %2d objective: %'lld", Iter, CurObjective);
    if (BestObjective == NOT_SET64 || BestObjective > CurObjective) {
      BestObjective = CurObjective;
      collectBestBuckets(DocumentBegin, DocumentEnd);
    }
    // Try to stop early
    if (IterObjective.size() >= MaxEqualAttempts) {
      bool AllEqual = true;
      for (size_t R = IterObjective.size() - MaxEqualAttempts; R < IterObjective.size(); R++) {
        if (IterObjective[R] != CurObjective) {
          AllEqual = false;
          break;
        }
      }
      if (AllEqual)
        break;
    }
    IterObjective.push_back(CurObjective);
  }

  if (BestObjective != NOT_SET64)
    assignBestBuckets(DocumentBegin, DocumentEnd);

  if (verbose >= 2 && Config.SplitAttempts > 1) {
    initializeSignatures(Signatures, DocumentBegin, DocumentEnd);
    LOG_IF(verbose >= 2, "    final   objective at bucket %d: %'lld", RootBucket, objective(Signatures, DocumentBegin, DocumentEnd));
  }
  Signatures.clear();

  //////////////////////////////////////////////////////////////////////////////
  // End of Optimization
  //////////////////////////////////////////////////////////////////////////////

  // Unset variables
  this->UseStrongMoveGains = false;
  this->LeftBucket = NOT_SET32;
  this->RightBucket = NOT_SET32;

  // Split documents wrt the resulting buckets
  const auto DocumentMid = std::stable_partition(
      DocumentBegin, DocumentEnd,
      [&](const Document *Doc) { return Doc->Bucket == LeftBucket; });

  const uint32_t MidOffset = Offset + std::distance(DocumentBegin, DocumentMid);
  const bool CanRecurse = Offset < MidOffset && MidOffset < Offset + NumDocuments;

  // Run recursively on the same thread
  if (CanRecurse) {
    CHECK(std::distance(DocumentBegin, DocumentMid) > 0);
    bisect(DocumentBegin, DocumentMid, RecDepth + 1, LeftBucket, Offset, true);

    CHECK(std::distance(DocumentMid, DocumentEnd) > 0);
    bisect(DocumentMid, DocumentEnd, RecDepth + 1, RightBucket, MidOffset, true);
  } else {
    order(DocumentBegin, DocumentEnd, Offset, Config.OrderAlg);
  }

  // Restore original utilities
  if (UseStrongMoveGains || Config.Backtrack) {
    restoreOrigUtilities(DocumentBegin, DocumentEnd);
    updateUtilities(DocumentBegin, DocumentEnd);
  }

  // Merge the two intervals optimally
  if (UseStrongMoveGains) {
    CHECK(isInitializedLB(DocumentBegin, DocumentEnd));
    mergeTwoIntervals(DocumentBegin, DocumentMid, DocumentEnd);

    // It is guaranteed that the document range consists of a single SCC
    if (Config.PostTuneInt) {
      tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O3, 0);
    } else {
      tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O1, 0);
    }
  }

  // Backtrack if the original solution is better
  if (Config.Backtrack) {
    CHECK(bucketsCorrect(DocumentBegin, DocumentEnd, true));
    uint64_t NewNumCrossings = countCrossings(DocumentBegin, DocumentEnd, NumUtilities);
    LOG_IF(verbose >= 2, "  completed the order at depth %d, bucket %d, new_crossings = %'lld", 
        RecDepth, RootBucket, NewNumCrossings);

    if (OrgNumCrossings < NewNumCrossings) {
      LOG_IF(verbose >= 2, "  reverting the order at depth %d with %d documents [OrgCrossings = %'lld; NewCrossings = %'lld]",
          RecDepth, NumDocuments, OrgNumCrossings, NewNumCrossings);
      if (UseStrongMoveGains) {
        CHECK(OrgIdOrder.size() == NumDocuments);
        std::unordered_map<uint32_t, uint32_t> Id2SortIndex;
        Id2SortIndex.reserve(NumDocuments);
        for (uint32_t I = 0; I < NumDocuments; I++) {
          CHECK(!Id2SortIndex.count(OrgIdOrder[I]));
          Id2SortIndex[OrgIdOrder[I]] = I;
        }
        for (uint32_t I = 0; I < NumDocuments; I++) {
          Document* Doc = DocumentBegin[I];
          CHECK(Id2SortIndex.count(Doc->id()));
          Doc->SortIdx = Id2SortIndex[Doc->id()];
        }
        std::sort(DocumentBegin, DocumentEnd, Document::SortIdxOrder());
        for (uint32_t I = 0; I < NumDocuments; I++) {
          Document* Doc = DocumentBegin[I];
          Doc->Bucket = Offset + I;
        }
        CHECK(bucketsCorrect(DocumentBegin, DocumentEnd, true));
      } else {
        order(DocumentBegin, DocumentEnd, Offset, OrgResult.second);
      }
    } else {
      LOG_IF(verbose >= 2, "  keeping the order at depth %d with %d documents [OrgCrossings = %'lld; NewCrossings = %'lld]",
          RecDepth, NumDocuments, OrgNumCrossings, NewNumCrossings);
    }
  }
}

template <class T>
void BalancedPartitioning<T>::runIterations(
    SignaturesType &Signatures,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd, 
    const uint32_t RecDepth) {
  // Initialize LB for move gains
  if (UseStrongMoveGains) {
    LOG_IF(verbose >= 3, "computeStrongMoveGains for |D| = %d", std::distance(DocumentBegin, DocumentEnd));
    computeStrongMoveGains(DocumentBegin, DocumentEnd);
  }

  // Run iterations
  uint32_t Iter = 0;
  while (Iter < Config.IterationsPerSplit) {
    LOG_IF(verbose >= 3, "    running iteration %d for depth %d", Iter, RecDepth);

    // Prepare move gains, if needed
    if (!UseStrongMoveGains)
      computeFastMoveGains(Signatures, DocumentBegin, DocumentEnd);

    uint32_t NumMovedDocuments = runIteration(Signatures, DocumentBegin, DocumentEnd, Iter);
    LOG_IF(verbose >= 3, "      moved %d documents", NumMovedDocuments);
    if (NumMovedDocuments == 0) {
      break;
    }
    Iter++;
  }
}

template <class T>
void BalancedPartitioning<T>::restoreOrigUtilities(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document *Doc = *It;
    Doc->restoreOrigEdges();
  }
}

uint32_t updateUtilities(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    bool UpdateOrigEdges) {
  CHECK(DocumentBegin != DocumentEnd);

  // Get the maximum utility adjacent to the given set of documents
  uint32_t MaxUtility = NOT_SET32;
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    const Document* Doc = *It;
    const uint32_t MaxVec = Doc->getEdges().back().utility;
    if (MaxUtility == NOT_SET32 || MaxUtility < MaxVec)
      MaxUtility = MaxVec;
  }
  CHECK(MaxUtility != NOT_SET32);
  const uint32_t NumUtilities = MaxUtility + 1;

  // Check which utilities are alive
  std::vector<bool> ActiveUtility(NumUtilities, false);
  uint32_t NumActive = 0;
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document *Doc = *It;
    for (const auto [U, _] : Doc->getEdges()) {
      if (!ActiveUtility[U]) {
        ActiveUtility[U] = true;
        NumActive++;
      }
    }
  }

  // Stop early if all are active
  if (NumActive == NumUtilities) {
    return NumUtilities;
  }

  // Count the (local) degree of each utility and compute their new indices
  std::vector<uint32_t> UtilityNodeIndex(NumUtilities);
  uint32_t CurIndex = 0;
  for (size_t I = 0; I < ActiveUtility.size(); I++) {
    if (ActiveUtility[I])
      UtilityNodeIndex[I] = CurIndex++;
  }

  // Update document adjacency lists
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document *Doc = *It;
    std::vector<EdgeTy> NewUtilityNodes;
    NewUtilityNodes.reserve(Doc->degree());
    for (const auto [U, W] : Doc->getEdges()) {
      CHECK(U <  NumUtilities && ActiveUtility[U]);
      CHECK(UtilityNodeIndex[U] < NumActive);
      NewUtilityNodes.emplace_back(UtilityNodeIndex[U], W);
    }
    if (UpdateOrigEdges)
      Doc->setOrigEdges(NewUtilityNodes);
    else
      Doc->setEdges(NewUtilityNodes);
  }

  // LOG_IF(UpdateOrigEdges, "removed %d unusued utilities", NumUtilities - NumActive);
  return NumActive;
}

template <class T>
void BalancedPartitioning<T>::initializeSignatures(
    SignaturesType &Signatures,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  for (auto It = Signatures.begin(); It != Signatures.end(); It++) {
    UtilitySignature& Signature = *It;
    Signature.LeftCount = 0;
    Signature.RightCount = 0;
  }

  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document *Doc = *It;
    // Initialize LastMovedIter
    Doc->LastMovedIter = NOT_SET32;
    if (Doc->Bucket == LeftBucket) {
      for (const auto [U, W] : Doc->getEdges()) {
        UtilitySignature& Signature = Signatures[U];
        Signature.LeftCount += W;
      }
    } else {
      for (const auto [U, W] : Doc->getEdges()) {
        UtilitySignature& Signature = Signatures[U];
        Signature.RightCount += W;
      }
    }
  }
}

template <class T>
void BalancedPartitioning<T>::initializeSignatureFastGains(
    SignaturesType &Signatures,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  // Compute Y->LB and X->RB
  int64_t XRB = 0;
  int64_t XLB = 0;
  for (auto It = Signatures.begin(); It != Signatures.end(); It++) {
    UtilitySignature& Signature = *It;

    Signature.XRB = XRB;
    Signature.XLB = XLB;
    XRB += Signature.RightCount;
    XLB += Signature.LeftCount;
  }
  int64_t YLB = 0;
  int64_t YRB = 0;
  for (auto It = Signatures.rbegin(); It != Signatures.rend(); It++) {
    UtilitySignature& Signature = *It;

    Signature.YLB = YLB;
    Signature.YRB = YRB;
    YLB += Signature.LeftCount;
    YRB += Signature.RightCount;
  }
}

template <class T>
void BalancedPartitioning<T>::initializeLBMatrix(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  CHECK(UseStrongMoveGains);
  CHECK(0 < NumDocuments && NumDocuments <= Config.MaxDocsLB);
  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));

  for (uint32_t I = 0; I < NumDocuments; I++) {
    Document* Doc = DocumentBegin[I];
    Doc->AdjLB.clear();
    Doc->AdjLB.reserve(NumDocuments);
  }

  // Make sure LB values are pre-computed
  for (uint32_t L = 0; L < NumDocuments; L++) {
    for (uint32_t R = L + 1; R < NumDocuments; R++) {
      Document* DocL = DocumentBegin[L];
      Document* DocR = DocumentBegin[R];
      const uint64_t LB_LR = LBData[DocumentBegin[L]->LocalIndex][DocumentBegin[R]->LocalIndex];
      const uint64_t LB_RL = LBData[DocumentBegin[R]->LocalIndex][DocumentBegin[L]->LocalIndex];

      if (LB_LR != LB_RL) {
        DocL->AdjLB.push_back(R);
        DocR->AdjLB.push_back(L);
      }
    }
  }
}

template <class T>
void BalancedPartitioning<T>::computeFastMoveGains(
    SignaturesType &Signatures,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  initializeSignatureFastGains(Signatures, DocumentBegin, DocumentEnd);

  for (UtilitySignature& Signature : Signatures) {
    CHECK(Signature.LeftCount > 0 || Signature.RightCount > 0);

    // Revise the objective?
    int64_t DiffCross = (Signature.XRB - Signature.YLB);
    int64_t DiffLocal = (Signature.XLB - Signature.YRB);

    Signature.CachedCostLR = DiffCross + DiffLocal;
  }

  // Move gains
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document *Doc = *It;

    int64_t Gain = 0;
    for (const auto [U, W] : Doc->getEdges()) {
      const UtilitySignature& Signature = Signatures[U];
      Gain += Signature.CachedCostLR * int64_t(W);
    }
    if (Doc->Bucket != LeftBucket) {
      Gain = -Gain;
    }

    Doc->MoveGain = Gain;
  }
}

int64_t calcObjective(const uint64_t LB_LR, const uint64_t LB_RL) {
  constexpr int64_t Sum = 100;
  constexpr int64_t MinScale = 95;
  constexpr int64_t MaxScale = Sum - MinScale;
  const int64_t MinLB = int64_t(std::min(LB_LR, LB_RL));
  const int64_t MaxLB = int64_t(std::max(LB_LR, LB_RL));
  return int64_t(MinLB * MinScale + MaxLB * MaxScale) / Sum;
}

template <class T>
int64_t BalancedPartitioning<T>::strongMoveGain(
    uint64_t LB_LR, uint64_t LB_RL,
    const uint32_t Doc1Bucket, const uint32_t Doc2Bucket) const {
  constexpr uint64_t ForcedWeight = 1000;
  if (LB_LR == 0 && LB_RL > 0)
    LB_RL *= ForcedWeight;
  else if (LB_RL == 0 && LB_LR > 0)
    LB_LR *= ForcedWeight;

  const int64_t MinObj = calcObjective(LB_LR, LB_RL);

  if (Doc1Bucket == Doc2Bucket) {
    return MinObj - int64_t(Doc1Bucket == LeftBucket ? LB_RL : LB_LR);
  }
  return int64_t(Doc1Bucket == LeftBucket ? LB_LR : LB_RL) - MinObj;
}

template <class T>
void BalancedPartitioning<T>::computeStrongMoveGains(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  CHECK(UseStrongMoveGains);

  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  // Compute move gains for every document
  for (uint32_t L = 0; L < NumDocuments; L++) {
    Document* Doc1 = DocumentBegin[L];
    Doc1->MoveGain = 0;

    for (uint32_t R : Doc1->AdjLB) {
      CHECK(L != R);
      const Document* Doc2 = DocumentBegin[R];

      const uint64_t LB_LR = LBData[Doc1->LocalIndex][Doc2->LocalIndex];
      const uint64_t LB_RL = LBData[Doc2->LocalIndex][Doc1->LocalIndex];

      CHECK(LB_LR != LB_RL);

      // Get the move gain
      Doc1->MoveGain += strongMoveGain(LB_LR, LB_RL, Doc1->Bucket, Doc2->Bucket);
    }
  }
}

template <class T>
void BalancedPartitioning<T>::updateStrongMoveGains(
    Document *Doc, uint32_t DocIndex,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
  CHECK(UseStrongMoveGains);

  // Compute move gains for every document
  const uint32_t CurBucket = Doc->Bucket;
  const uint32_t PrvBucket = Doc->Bucket == LeftBucket ? RightBucket : LeftBucket;

  for (uint32_t R : Doc->AdjLB) {
    Document* Doc2 = DocumentBegin[R];

    const uint64_t LB_LR = LBData[Doc->LocalIndex][Doc2->LocalIndex];
    const uint64_t LB_RL = LBData[Doc2->LocalIndex][Doc->LocalIndex];

    // Update for Doc1
    Doc->MoveGain -= strongMoveGain(LB_LR, LB_RL, PrvBucket, Doc2->Bucket);
    Doc->MoveGain += strongMoveGain(LB_LR, LB_RL, CurBucket, Doc2->Bucket);

    // Update for Doc2
    Doc2->MoveGain -= strongMoveGain(LB_RL, LB_LR, Doc2->Bucket, PrvBucket);
    Doc2->MoveGain += strongMoveGain(LB_RL, LB_LR, Doc2->Bucket, CurBucket);
  }
}

template <class T>
uint32_t BalancedPartitioning<T>::runIteration(
    SignaturesType &Signatures,
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd, 
    const uint32_t Iter) {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);

  // Compute move gains
  typedef std::pair<int64_t, uint32_t> GainPair;
  std::vector<GainPair> LeftGains;
  LeftGains.reserve(NumDocuments);
  std::vector<GainPair> RightGains;
  RightGains.reserve(NumDocuments);
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    const Document *Doc = *It;
    const uint32_t Index = It - DocumentBegin;
    const int64_t Gain = Doc->MoveGain;
    if (Doc->Bucket == LeftBucket)
      LeftGains.push_back({Gain, Index});
    else
      RightGains.push_back({Gain, Index});
  }

  // Sort gains
  std::sort(LeftGains.begin(), LeftGains.end(), std::greater<GainPair>());
  std::sort(RightGains.begin(), RightGains.end(), std::greater<GainPair>());

  // Exchange buckets and update signatures
  uint32_t NumMovedDocs = 0;
  uint32_t MinSize = std::min(LeftGains.size(), RightGains.size());
  for (uint32_t I = 0; I < MinSize; I++) {
    if (LeftGains[I].first + RightGains[I].first <= 0)
      break;

    Document* DocL = DocumentBegin[LeftGains[I].second];
    Document* DocR = DocumentBegin[RightGains[I].second];

    // Try to swap the two documents
    const bool UseApproxGainL = (!Config.UseActualGains && AbsDiff(DocL->LastMovedIter, Iter) > 1 && Iter < 15);
    const int64_t DocLGain = UseApproxGainL ? LeftGains[I].first : DocL->MoveGain;
    if (DocLGain > 0 && moveDocument(DocL, Signatures)) {
      NumMovedDocs += 1;
      DocL->LastMovedIter = Iter;
      if (UseStrongMoveGains)
        updateStrongMoveGains(DocL, LeftGains[I].second, DocumentBegin, DocumentEnd);
    }

    const bool UseApproxGainR = (!Config.UseActualGains && AbsDiff(DocR->LastMovedIter, Iter) > 1 && Iter < 15);
    const int64_t DocRGain = UseApproxGainR ? RightGains[I].first : DocR->MoveGain;
    if (DocRGain > 0 && moveDocument(DocR, Signatures)) {
      NumMovedDocs += 1;
      DocR->LastMovedIter = Iter;
      if (UseStrongMoveGains)
        updateStrongMoveGains(DocR, RightGains[I].second, DocumentBegin, DocumentEnd);
    }
  }

  return NumMovedDocs;
}

template <class T>
bool BalancedPartitioning<T>::moveDocument(Document *Doc,
                                        SignaturesType &Signatures) const {
  // Sometimes we skip the move. This helps to escape local optima
  if (Rand::nextDouble(RNG) <= Config.SkipProbability)
    return false;

  // Update the current bucket and all signatures
  if (Doc->Bucket == LeftBucket) {
    Doc->Bucket = RightBucket;
    for (const auto [U, W] : Doc->getEdges()) {
      UtilitySignature& Signature = Signatures[U];
      Signature.LeftCount -= W;
      Signature.RightCount += W;
    }
  } else {
    Doc->Bucket = LeftBucket;
    for (const auto [U, W] : Doc->getEdges()) {
      UtilitySignature& Signature = Signatures[U];
      Signature.LeftCount += W;
      Signature.RightCount -= W;
    }
  }

  return true;
}

template <class T>
void BalancedPartitioning<T>::order(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    uint32_t StartBucket, OrderAlgT OrderAlgorithm) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  if (NumDocuments <= Config.LeafInterval)
    LOG_IF(verbose >= 2 && NumDocuments > 1, "    ordering %d docs with opt-%d", NumDocuments, Config.LeafInterval);
  else
    LOG_IF(verbose >= 2 && NumDocuments > 1, "    ordering %d docs with alg-%d", NumDocuments, (int)OrderAlgorithm);

  auto assignBuckets = [&](uint32_t MinBucket) {
    for (auto It = DocumentBegin; It != DocumentEnd; It++) {
      (*It)->Bucket = MinBucket++;
    }
  };

  auto computeLocal = [&]() {
    std::vector<uint32_t> UNodes;
    for (auto It = DocumentBegin; It != DocumentEnd; It++) {
      Document* Doc = *It;
      UNodes.clear();
      UNodes.reserve(Doc->degree());
      for (const auto [U, _] : Doc->getEdges()) {
        UNodes.push_back(U);
      }
      Doc->LocalMedian = median(UNodes);
      Doc->LocalAverage = average(UNodes);
    }
  };

  if (OrderAlgorithm == OrderAlgT::NONE || NumDocuments == 1) {
    // pass
  } else if (OrderAlgorithm == OrderAlgT::INPUT) {
    std::sort(DocumentBegin, DocumentEnd, Document::IdOrder());
  } else if (NumDocuments <= Config.LeafInterval) {
    std::sort(DocumentBegin, DocumentEnd, Document::MedianOrder());
    assignBuckets(StartBucket);
    initLB(DocumentBegin, DocumentEnd);
    sortIntervalOpt(DocumentBegin, DocumentEnd, false);
  } else if (OrderAlgorithm == OrderAlgT::RANDOM) {
    std::vector<uint32_t> RndPerm = Rand::permutation(NumDocuments, RNG);
    for (uint32_t I = 0; I < NumDocuments; I++) {
      Document* Doc = DocumentBegin[I];
      Doc->SortIdx = RndPerm[I];
    }
    std::sort(DocumentBegin, DocumentEnd, Document::SortIdxOrder());
  } else if (OrderAlgorithm == OrderAlgT::MEDIAN) {
    std::sort(DocumentBegin, DocumentEnd, Document::MedianOrder());
  } else if (OrderAlgorithm == OrderAlgT::AVERAGE) {
    std::sort(DocumentBegin, DocumentEnd, Document::AverageOrder());
  } else if (OrderAlgorithm == OrderAlgT::LOCAL_MEDIAN || OrderAlgorithm == OrderAlgT::LOCAL_AVERAGE) {
    computeLocal();
    if (OrderAlgorithm == OrderAlgT::LOCAL_MEDIAN) {
      std::sort(DocumentBegin, DocumentEnd, Document::LocalMedianOrder());
    } else {
      std::sort(DocumentBegin, DocumentEnd, Document::LocalAverageOrder());
    }
  } else if (OrderAlgorithm == OrderAlgT::AVERAGE_OPT) {
    CHECK(NumDocuments <= Config.MaxDocsLB);
    std::sort(DocumentBegin, DocumentEnd, Document::AverageOrder());
    assignBuckets(StartBucket);
    tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O3, 0);
  } else if (OrderAlgorithm == OrderAlgT::MEDIAN_OPT) {
    CHECK(NumDocuments <= Config.MaxDocsLB);
    computeLocal();
    std::sort(DocumentBegin, DocumentEnd, Document::LocalMedianOrder());
    assignBuckets(StartBucket);
    tuneInterval(DocumentBegin, DocumentEnd, false, OptLevelT::O3, 0);
  } else
    ERROR("unknown OrderAlg");

  // Assign buckets
  assignBuckets(StartBucket);
}

template <class T>
std::pair<uint64_t, OrderAlgT>
BalancedPartitioning<T>::orderForBacktrack(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    const uint32_t StartBucket, const uint32_t NumUtilities,
    std::vector<uint32_t>& OrgIdOrder) const {
  auto saveBuckets = [&]() {
    if (!UseStrongMoveGains)
      return;
    OrgIdOrder.clear();
    OrgIdOrder.reserve(std::distance(DocumentBegin, DocumentEnd));
    for (auto It = DocumentBegin; It != DocumentEnd; It++) {
      OrgIdOrder.push_back((*It)->id());
    }
  };

  uint64_t BestCrossings = MAX_64;
  OrderAlgT BestAlg;

  // try LOCAL_MEDIAN
  order(DocumentBegin, DocumentEnd, StartBucket, OrderAlgT::LOCAL_MEDIAN);
  const uint64_t numCrossLM = countCrossings(DocumentBegin, DocumentEnd, NumUtilities);
  if (numCrossLM < BestCrossings) {
    BestCrossings = numCrossLM;
    BestAlg = OrderAlgT::LOCAL_MEDIAN;
    saveBuckets();
  }
  if (numCrossLM == 0) {
    return {numCrossLM, OrderAlgT::LOCAL_MEDIAN};
  }

  // try LOCAL_AVERAGE
  order(DocumentBegin, DocumentEnd, StartBucket, OrderAlgT::LOCAL_AVERAGE);
  const uint64_t numCrossLA = countCrossings(DocumentBegin, DocumentEnd, NumUtilities);
  if (numCrossLA < BestCrossings) {
    BestCrossings = numCrossLA;
    BestAlg = OrderAlgT::LOCAL_AVERAGE;
    saveBuckets();
  }

  const bool UseHeavy = Rand::check(0.1, RNG);

  // try MEDIAN_OPT
  if (UseStrongMoveGains && UseHeavy) {
    order(DocumentBegin, DocumentEnd, StartBucket, OrderAlgT::MEDIAN_OPT);
    const uint64_t numCrossMO = countCrossings(DocumentBegin, DocumentEnd, NumUtilities);
    if (numCrossMO < BestCrossings) {
      BestCrossings = numCrossMO;
      BestAlg = OrderAlgT::MEDIAN_OPT;
      saveBuckets();
    }
  }

  // try AVERAGE_OPT
  if (UseStrongMoveGains && UseHeavy) {
    order(DocumentBegin, DocumentEnd, StartBucket, OrderAlgT::AVERAGE_OPT);
    const uint64_t numCrossAO = countCrossings(DocumentBegin, DocumentEnd, NumUtilities);
    if (numCrossAO < BestCrossings) {
      BestCrossings = numCrossAO;
      BestAlg = OrderAlgT::AVERAGE_OPT;
      saveBuckets();
    }
  }

  CHECK(BestCrossings != MAX_64);
  return {BestCrossings, BestAlg};
}

template <class T>
void BalancedPartitioning<T>::split(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    const uint32_t StartBucket) const {
  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  CHECK(NumDocuments > 0 && "incorrect number of documents");

  uint32_t DocsPerBucket = (NumDocuments + 1) / 2;
  CHECK(DocsPerBucket > 0 && "incorrect number of docs per bucket");

  auto assignBuckets = [&]() {
    auto It = DocumentBegin;
    for (uint32_t I = 0; I < 2; I++) {
      for (uint32_t J = 0; J < DocsPerBucket && It != DocumentEnd; J++) {
        (*It)->Bucket = StartBucket + I;
        It++;
      }
    }
  };

  if (Config.SplitAlg == SplitAlgT::INPUT) {
    CHECK(Config.SplitAttempts == 1, "multi-splits can be applied only for SplitAlgT::RANDOM");
    std::nth_element(DocumentBegin, DocumentBegin + DocsPerBucket, DocumentEnd,
                    [](const Document *DocL, const Document *DocR) {
                      return DocL->id() < DocR->id();
                    });
    assignBuckets();
  } else if (Config.SplitAlg == SplitAlgT::MEDIAN) {
    CHECK(Config.SplitAttempts == 1, "multi-splits can be applied only for SplitAlgT::RANDOM");
    std::nth_element(DocumentBegin, DocumentBegin + DocsPerBucket, DocumentEnd,
                    [](const Document *DocL, const Document *DocR) {
                      return std::tuple(DocL->AllMedian, DocL->AllAverage, DocL->id()) <
                             std::tuple(DocR->AllMedian, DocR->AllAverage, DocR->id());
                    });
    assignBuckets();
  } else if (Config.SplitAlg == SplitAlgT::AVERAGE) {
    CHECK(Config.SplitAttempts == 1, "multi-splits can be applied only for SplitAlgT::RANDOM");
    std::nth_element(DocumentBegin, DocumentBegin + DocsPerBucket, DocumentEnd,
                    [](const Document *DocL, const Document *DocR) {
                      return std::tuple(DocL->AllAverage, DocL->AllMedian, DocL->id()) <
                             std::tuple(DocR->AllAverage, DocR->AllMedian, DocR->id());
                    });
    assignBuckets();
  } else if (Config.SplitAlg == SplitAlgT::RANDOM) {
    std::vector<uint32_t> RndPerm = Rand::permutation(NumDocuments, RNG);
    uint32_t CurIdx = 0;
    for (uint32_t I = 0; I < 2; I++) {
      for (uint32_t J = 0; J < DocsPerBucket && CurIdx < NumDocuments; J++) {
        Document* Doc = DocumentBegin[RndPerm[CurIdx]];
        Doc->Bucket = StartBucket + I;
        CurIdx++;
      }
    }
  } else {
    ERROR("unknown SplitAlg");
  }
}

template <class T>
uint64_t BalancedPartitioning<T>::countCrossings(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t NumUtilities) const {
  const uint32_t n = NumUtilities;
  CHECK(stree.size() >= 2 * n);
  std::fill(stree.begin(), stree.begin() + 2 * n, 0);

  uint64_t NumCrossings = 0;
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    const Document* Doc = *It;
    // Count crossings
    for (const auto [U, W] : Doc->getEdges()) {
      NumCrossings += sumST(n, stree, U + 1, n) * W;
    }
    // Self-crossings: This might return an incorrect result in the presence of
    // multi-edges. Ignore for now
    NumCrossings += Doc->getSelfCrossings();
    // Update the tree
    for (const auto [U, W] : Doc->getEdges()) {
      addST(n, stree, U, W);
    }
  }
  return NumCrossings;
}

template <class T>
uint64_t BalancedPartitioning<T>::countCrossings(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  return countCrossings(DocumentBegin, DocumentEnd, MaxUtility + 1);
}

template <class T>
uint64_t BalancedPartitioning<T>::objective(
    SignaturesType &Signatures,
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  if (!UseStrongMoveGains) {
    initializeSignatureFastGains(Signatures, DocumentBegin, DocumentEnd);

    uint64_t ExpCrossings = 0;
    for (auto It = DocumentBegin; It != DocumentEnd; It++) {
      const Document* Doc = *It;
      if (Doc->Bucket == LeftBucket) {
        for (const auto [U, W] : Doc->getEdges()) {
          const UtilitySignature& Signature = Signatures[U];
          ExpCrossings += uint64_t(Signature.XLB + Signature.YLB + 2 * Signature.XRB) * W;
        }
      } else {
        for (const auto [U, W] : Doc->getEdges()) {
          const UtilitySignature& Signature = Signatures[U];
          ExpCrossings += uint64_t(Signature.XRB + Signature.YRB + 2 * Signature.YLB) * W;
        }
      }
    }
    return ExpCrossings;
  }

  CHECK(UseStrongMoveGains);
  uint64_t NumCrossings = 0;
  for (auto It1 = DocumentBegin; It1 != DocumentEnd; It1++) {
    const Document* Doc1 = *It1;
    for (auto It2 = It1 + 1; It2 != DocumentEnd; It2++) {
      const Document* Doc2 = *It2;
      if (Doc1->Bucket == Doc2->Bucket) {
        const uint64_t LB_LR = LBData[Doc1->LocalIndex][Doc2->LocalIndex];
        const uint64_t LB_RL = LBData[Doc2->LocalIndex][Doc1->LocalIndex];
        // NumCrossings += std::min(LB_LR, LB_RL);
        NumCrossings += calcObjective(LB_LR, LB_RL);
      } else if (Doc1->Bucket < Doc2->Bucket) {
        NumCrossings += LBData[Doc1->LocalIndex][Doc2->LocalIndex];
      } else {
        NumCrossings += LBData[Doc2->LocalIndex][Doc1->LocalIndex];
      }
    }
  }
  return NumCrossings;
}

template <class T>
void BalancedPartitioning<T>::collectBestBuckets(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document* Doc = *It;
    Doc->BestBucket = Doc->Bucket;
  }
}

template <class T>
void BalancedPartitioning<T>::assignBestBuckets(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd) const {
  for (auto It = DocumentBegin; It != DocumentEnd; It++) {
    Document* Doc = *It;
    Doc->Bucket = Doc->BestBucket;
  }
}

template <class T>
void BalancedPartitioning<T>::sortByOrder(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const std::vector<uint32_t>& Order) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N == Order.size() && N > 0);
#ifdef DEBUG
  std::vector<uint32_t> inv = inverse(Order);
  CHECK(inv.size() == N);
#endif

  uint32_t MinBucket = DocumentBegin[0]->Bucket;
  uint32_t MaxBucket = DocumentBegin[0]->Bucket;
  for (uint32_t I = 0; I < N; I++) {
    const Document* Doc = DocumentBegin[I];
    MinBucket = std::min(MinBucket, Doc->Bucket);
    MaxBucket = std::max(MaxBucket, Doc->Bucket);
  }
  if (MinBucket == MaxBucket) {
    for (uint32_t I = 0; I < N; I++) {
      Document* Doc = *(DocumentBegin + Order[I]);
      Doc->SortIdx = I;
    }
    std::sort(DocumentBegin, DocumentEnd,
              [](const Document *DocL, const Document *DocR) {
                return DocL->SortIdx < DocR->SortIdx;
              });
    return;
  }

  if (MaxBucket - MinBucket + 1 != N) {
    LOG_IF(verbose, "N = %d", N);
    for (uint32_t I = 0; I < N; I++) {
      const Document* Doc = DocumentBegin[I];
      LOG_IF(verbose, "  %d", Doc->Bucket);
    }
  }

  CHECK(MaxBucket - MinBucket + 1 == N, "MinBucket = %d, MaxBucket = %d, N = %d", MinBucket, MaxBucket, N);
  for (uint32_t I = 0; I < N; I++) {
    Document* Doc = *(DocumentBegin + Order[I]);
    CHECK(MinBucket <= Doc->Bucket && Doc->Bucket <= MaxBucket);
    Doc->Bucket = MinBucket + I;
  }
  std::sort(DocumentBegin, DocumentEnd, Document::BucketOrder());
}

template <class T>
uint64_t BalancedPartitioning<T>::computeLowerBound(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const uint32_t IntervalSize) const {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  CHECK(N <= Config.MaxDocsLB);
  CHECK(IntervalSize >= 2);
  initLB(DocumentBegin, DocumentEnd);

  std::vector<uint64_t> OptLB(N, 0);
  std::vector<uint32_t> Parent(N, NOT_SET32);

  auto getIntervalLB = [&](uint32_t L, uint32_t R) {
    CHECK(L < R);
    CHECK(R - L + 1 <= IntervalSize);
    const std::vector<Document*>::iterator DocBegin = DocumentBegin + L;
    const std::vector<Document*>::iterator DocEnd = DocumentBegin + R + 1;

    const auto [CurCrossings, OptCrossings] = computeCurOptCrossings(DocBegin, DocEnd);
    CHECK(CurCrossings >= OptCrossings);
    return CurCrossings - OptCrossings;
  };

  for (uint32_t I = 1; I < N; I++) {
    const uint32_t Start = I + 1 >= IntervalSize ? I + 1 - IntervalSize : 0;
    for (uint32_t J = Start; J < I; J++) {
      const uint64_t Res = OptLB[J] + getIntervalLB(J, I);
      if (Parent[I] == NOT_SET32 || OptLB[I] < Res) {
        OptLB[I] = Res;
        Parent[I] = J;
      }
    }
  }

  std::vector<std::pair<uint32_t, uint32_t>> TakenIntervals;
  uint32_t Now = N - 1;
  while (Parent[Now] != NOT_SET32) {
    uint32_t Prev = Parent[Now];
    CHECK(Prev < Now);
    TakenIntervals.push_back({Prev, Now});
    Now = Prev;
  }
  CHECK(Now == 0);

  std::vector<std::vector<bool>> Taken(N, std::vector<bool>(N));
  uint64_t IntervalLB = 0;
  for (const auto& [Begin, End] : TakenIntervals) {
    for (uint32_t L = Begin; L <= End; L++) {
      for (uint32_t R = L + 1; R <= End; R++) {
        CHECK(!Taken[L][R] && !Taken[R][L]);
        const uint64_t LB_LR = LBData[DocumentBegin[L]->LocalIndex][DocumentBegin[R]->LocalIndex];

        IntervalLB += LB_LR;
        Taken[L][R] = Taken[R][L] = true;
      }
    }
  }

  for (uint32_t L = 0; L < N; L++) {
    for (uint32_t R = L + 1; R < N; R++) {
      if (Taken[L][R])
        continue;
      const uint64_t LB_LR = LBData[DocumentBegin[L]->LocalIndex][DocumentBegin[R]->LocalIndex];
      const uint64_t LB_RL = LBData[DocumentBegin[R]->LocalIndex][DocumentBegin[L]->LocalIndex];

      IntervalLB += std::min(LB_LR, LB_RL);
    }
  }

  // Account for (merged) document duplicates
  for (uint32_t L = 0; L < N; L++) {
    const Document* Doc = DocumentBegin[L];
    IntervalLB += Doc->getSelfCrossings();
  }

  // Sanity check: count the actual number of crossings (this isn't guaranteed in certain cases)
  // uint64_t ActualCrossings = 0;
  // for (uint32_t L = 0; L < N; L++) {
  //   for (uint32_t R = L + 1; R < N; R++) {
  //     ActualCrossings += LBM.get(DocumentBegin[L], DocumentBegin[R]);
  //   }
  // }
  // for (uint32_t L = 0; L < N; L++) {
  //   const Document* Doc = DocumentBegin[L];
  //   ActualCrossings += Doc->getSelfCrossings();
  // }
  // CHECK(ActualCrossings >= IntervalLB);

  return IntervalLB;
}

template <class T>
void BalancedPartitioning<T>::computeLowerBound(
    std::vector<Document *>& Documents,
    uint32_t IntervalSize) {
  CHECK(bucketsCorrect(Documents.begin(), Documents.end()));
  CHECK(Documents.size() <= Config.MaxDocsLB);

  uint64_t IntervalLB = computeLowerBound(Documents.begin(), Documents.end(), IntervalSize);

  LOG_IF(verbose >= 1, " interval-%d lower bound: %'lld", IntervalSize, IntervalLB);
  if (LowerBound == NOT_SET64 || LowerBound < IntervalLB) {
    LowerBound = IntervalLB;
  }
}

bool bucketsCorrect(
    const std::vector<Document*>::iterator& DocumentBegin,
    const std::vector<Document*>::iterator& DocumentEnd,
    const bool StrictOrder) {
#ifdef DEBUG
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  if (N == 0)
    return true;

  const uint32_t MinBucket = DocumentBegin[0]->Bucket;
  const uint32_t MaxBucket = DocumentBegin[N - 1]->Bucket;
  if (MinBucket == MaxBucket && !StrictOrder)
    return true;

  if (MaxBucket - MinBucket + 1 != N) {
    //LOG_IF(verbose, "1: N = %d, MinBucket = %d, MaxBucket = %d", N, MinBucket, MaxBucket);
    // for (uint32_t I = 1; I < N; I++) {
    //   LOG_IF(verbose, "   Bucket[%d] = %d", I, DocumentBegin[I]->Bucket);
    // }
    return false;
  }
  for (uint32_t I = 1; I < N; I++) {
    if (DocumentBegin[I]->Bucket != MinBucket + I) {
      //LOG_IF(verbose, "2: N = %d, MinBucket = %d, MaxBucket = %d", N, MinBucket, MaxBucket);
      return false;
    }
  }
#endif
  return true;
}

std::vector<uint32_t> extractOrderedIds(const std::vector<Document *>& Docs) {
  // Verify that every document gets a correct bucket
  bool IsSorted = true;
  std::vector<uint32_t> DocOrder(Docs.size());
  for (size_t I = 0; I < Docs.size(); I++) {
    const Document* Doc = Docs[I];
    if (Doc->Bucket != I)
      IsSorted = false;
    DocOrder[I] = I;
    // CHECK(0 <= Doc->Bucket && Doc->Bucket < Docs.size());
  }
  if (!IsSorted) {
    // Sort documents by the resulting buckets
    std::sort(DocOrder.begin(), DocOrder.end(), [&](uint32_t L, const uint32_t R) {
      return Docs[L]->Bucket < Docs[R]->Bucket;
    });
  }
  std::vector<uint32_t> Order;
  Order.reserve(Docs.size());
  for (size_t I = 0; I < Docs.size(); I++) {
    // Verify that buckets are distinct
    CHECK(I + 1 == Docs.size() || Docs[DocOrder[I]]->Bucket < Docs[DocOrder[I + 1]]->Bucket);
    const Document *Doc = Docs[DocOrder[I]];
    for (uint32_t Id : Doc->ids()) {
      Order.push_back(Id);
    }
  }
  return Order;
}

void collectBuckets(
    const std::vector<Document *>& Documents,
    std::vector<std::pair<uint32_t, uint32_t>>& Buckets) {
  Buckets.clear();
  Buckets.reserve(Documents.size());
  for (uint32_t I = 0; I < Documents.size(); I++) {
    Buckets.push_back({Documents[I]->id(), Documents[I]->Bucket});
  }
  std::sort(Buckets.begin(), Buckets.end());
}

void assignBuckets(
    std::vector<Document *>& Documents,
    const std::vector<std::pair<uint32_t, uint32_t>>& Buckets) {
  CHECK(Documents.size() == Buckets.size());
  std::sort(Documents.begin(), Documents.end(), Document::IdOrder());
  for (uint32_t I = 0; I < Documents.size(); I++) {
    Document* Doc = Documents[I];
    CHECK(Doc->id() == Buckets[I].first);
    Doc->Bucket = Buckets[I].second;
  }
  std::sort(Documents.begin(), Documents.end(), Document::BucketOrder());
}

template class BalancedPartitioning<uint8_t>;
template class BalancedPartitioning<uint16_t>;
template class BalancedPartitioning<uint32_t>;
template class BalancedPartitioning<uint64_t>;
