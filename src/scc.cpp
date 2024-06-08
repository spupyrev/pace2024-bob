#include "bp.h"
#include "logging.h"
#include "common.h"

#include <queue>

void dfsForward(const std::vector<Document *>::iterator &DocumentBegin,
                const std::vector<Document *>::iterator &DocumentEnd,
                uint32_t Now,
                std::vector<bool>& Used,
                std::vector<uint32_t>& Order) {
  const Document* Doc = DocumentBegin[Now];
  Used[Now] = true;

  for (uint32_t Next : Doc->SuccLB) {
    if (Used[Next])
      continue;
    dfsForward(DocumentBegin, DocumentEnd, Next, Used, Order);
  }

  Order.push_back(Now);
}

void dfsBackward(const std::vector<Document *>::iterator &DocumentBegin,
                 const std::vector<Document *>::iterator &DocumentEnd,
                 uint32_t Now,
                 std::vector<bool>& Used,
                 std::vector<uint32_t>& Component) {
  const Document* Doc = DocumentBegin[Now];
  Used[Now] = true;
  Component.push_back(Now);

  for (uint32_t Next : Doc->PredLB) {
    if (Used[Next])
      continue;
    dfsBackward(DocumentBegin, DocumentEnd, Next, Used, Component);
  }
}

std::vector<std::vector<uint32_t>>
stronglyConnectedComponents(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) {
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);

  std::vector<bool> Used(N);
  std::vector<uint32_t> Order;
  Order.reserve(N);
  std::vector<uint32_t> Component;
  Component.reserve(N);

  Used.assign(N, false);

  for (uint32_t I = 0; I < N; I++) {
    if (Used[I])
      continue;
    dfsForward(DocumentBegin, DocumentEnd, I, Used, Order);
  }
  CHECK(Order.size() == N);
  std::reverse(Order.begin(), Order.end());

  Used.assign(N, false);
  std::vector<std::vector<uint32_t>> Result;
  for (uint32_t I : Order) {
    if (Used[I])
      continue;
    dfsBackward(DocumentBegin, DocumentEnd, I, Used, Component);
    Result.push_back(Component);
    Component.clear();
  }
  return Result;
}

template <class T>
uint64_t BalancedPartitioning<T>::splitAndReorderBySCCs(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd,
    std::vector<IterRange>& Result) const {
  Result.clear();
  const uint32_t N = std::distance(DocumentBegin, DocumentEnd);
  if (N > Config.MaxDocsLB || N <= 1) {
    Result.push_back({DocumentBegin, DocumentEnd});
    return 0;
  }

  CHECK(isInitializedLB(DocumentBegin, DocumentEnd));
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));

  // Don't run it if we already have the optimum
  const auto [CurCrossings, OptCrossings] = computeCurOptCrossings(DocumentBegin, DocumentEnd);
  if (CurCrossings == OptCrossings) {
    Result.push_back({DocumentBegin, DocumentEnd});
    return 0;
  }

  // Init reachability
  for (uint32_t I = 0; I < N; I++) {
    Document* Doc = DocumentBegin[I];
    Doc->SuccLB.clear();
    Doc->SuccLB.reserve(N);
    Doc->PredLB.clear();
    Doc->PredLB.reserve(N);
  }
  for (uint32_t L = 0; L < N; L++) {
    for (uint32_t R = L + 1; R < N; R++) {
      Document* DocL = DocumentBegin[L];
      Document* DocR = DocumentBegin[R];
      const uint64_t LB_LR = LBData[DocL->LocalIndex][DocR->LocalIndex];
      const uint64_t LB_RL = LBData[DocR->LocalIndex][DocL->LocalIndex];

      if (LB_LR < LB_RL) {
        DocL->SuccLB.push_back(R);
        DocR->PredLB.push_back(L);
      } else if (LB_LR > LB_RL) {
        DocL->PredLB.push_back(R);
        DocR->SuccLB.push_back(L);
      }
    }
  }

  // Find SCCs topologically ordered
  auto SCCs = stronglyConnectedComponents(DocumentBegin, DocumentEnd);
  if (SCCs.size() == 1) {
    CHECK(SCCs[0].size() == N);
    Result.push_back({DocumentBegin, DocumentEnd});
    return 0;
  }

  // Stable sorting the documents
  std::vector<uint32_t> TSIndex(N, NOT_SET32);
  for (uint32_t I = 0; I < SCCs.size(); I++) {
    for (uint32_t J = 0; J < SCCs[I].size(); J++) {
      TSIndex[SCCs[I][J]] = I;
    }
  }
  std::vector<uint32_t> Order = identity(N);
  std::stable_sort(Order.begin(), Order.end(), [&](const uint32_t L, const uint32_t R) {
    return std::make_tuple(TSIndex[L], L) < std::make_tuple(TSIndex[R], R);
  });
  sortByOrder(DocumentBegin, DocumentEnd, Order);
  CHECK(bucketsCorrect(DocumentBegin, DocumentEnd));

  // Prepare the result
  uint32_t TotalSize = 0;
  for (uint32_t I = 0; I < SCCs.size(); I++) {
    const std::vector<Document *>::iterator DocBegin = DocumentBegin + TotalSize;
    const std::vector<Document *>::iterator DocEnd = DocumentBegin + TotalSize + SCCs[I].size();
    Result.push_back({DocBegin, DocEnd});
    TotalSize += SCCs[I].size();
  }
  CHECK(TotalSize == N);

  const auto [NewCrossings, _] = computeCurOptCrossings(DocumentBegin, DocumentEnd);
  CHECK(NewCrossings <= CurCrossings);
  const uint64_t Delta = CurCrossings - NewCrossings;

  return Delta;
}

template class BalancedPartitioning<uint8_t>;
template class BalancedPartitioning<uint16_t>;
template class BalancedPartitioning<uint32_t>;
template class BalancedPartitioning<uint64_t>;
