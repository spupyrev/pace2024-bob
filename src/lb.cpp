#include "bp.h"
#include "logging.h"
#include "common.h"

uint64_t computeCrossings(const Document* DocL, const Document* DocR) {
  const std::vector<EdgeTy>& UtilsL = DocL->getOrigEdges();
  const std::vector<EdgeTy>& UtilsR = DocR->getOrigEdges();
  uint64_t Res = 0;
  if (UtilsL.back().utility <= UtilsR.front().utility) {
    // no crossings
  } else if (UtilsL.front().utility > UtilsR.back().utility) {
    // everything crosses
    Res = DocL->weightedDegree() * DocR->weightedDegree();
  } else {
    uint32_t J = 0;
    uint64_t NumRUtils = 0;
    for (uint32_t I = 0; I < UtilsL.size(); I++) {
      while (J < UtilsR.size() && UtilsR[J].utility < UtilsL[I].utility) {
        NumRUtils += UtilsR[J].weight;
        J++;
      }
      Res += NumRUtils * UtilsL[I].weight;
    }
  }
  return Res;
}

template <class T>
void BalancedPartitioning<T>::initializeLB(std::vector<Document *>& Documents) {
  const uint32_t N = Documents.size();
  for (uint32_t I = 0; I < N; I++) {
    CHECK(Documents[I]->Index == I);
  }

  const uint32_t MaxDocsLB = std::max(Config.MaxDocsLB, uint32_t(32));
  CHECK(N <= MaxDocsLB || !UseDenseLBData);

  CHECK(LBData.empty() && CachedDocs.empty());
  LBData = std::vector<std::vector<T>>(MaxDocsLB, std::vector<T>(MaxDocsLB, T(-1)));
  CachedDocs.reserve(MaxDocsLB);

  if (N > MaxDocsLB) {
    LOG_IF(verbose >= 1, "skipping initializeLB (%d bytes) due to too many documents (%d > %d)", sizeof(T), N, MaxDocsLB);
    return;
  }

  LOG_IF(verbose >= 1, "InitializeLB (%d bytes) for %d documents", sizeof(T), N);
  for (uint32_t L = 0; L < N; L++) {
    Documents[L]->LocalIndex = L;
    for (uint32_t R = L + 1; R < N; R++) {
      const Document* DocL = Documents[L];
      const Document* DocR = Documents[R];

      const uint64_t valueLR = computeCrossings(DocL, DocR);
      CHECK(valueLR < std::numeric_limits<T>::max(), "ResLR = %'lld", valueLR);
      const uint64_t valueRL = computeCrossings(DocR, DocL);
      CHECK(valueRL < std::numeric_limits<T>::max(), "ResRL = %'lld", valueRL);

      LBData[L][R] = valueLR;
      LBData[R][L] = valueRL;
    }
  }
}

template <class T>
bool BalancedPartitioning<T>::isInitializedLB(
    const std::vector<Document *>::iterator &DocumentBegin,
    const std::vector<Document *>::iterator &DocumentEnd) const {
#ifdef DEBUG
  if (UseDenseLBData)
    return true;

  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  for (uint32_t I = 0; I < NumDocuments; I++) {
    const Document* Doc = DocumentBegin[I];
    if (CachedDocs.find(Doc->Index) == CachedDocs.end()) {
      // LOG("not cached for doc idx = %d; |CachedDocs| = %d", Doc->Index, CachedDocs.size());
      return false;
    }
  }
#endif
  return true;
}

template <class T>
void BalancedPartitioning<T>::initLB(
    const std::vector<Document *>::iterator& DocumentBegin,
    const std::vector<Document *>::iterator& DocumentEnd) const {
  if (UseDenseLBData)
    return;

  const uint32_t NumDocuments = std::distance(DocumentBegin, DocumentEnd);
  LOG_IF(verbose >= 2, "initLB for %d docs [%d .. %d]",
      NumDocuments, DocumentBegin[0]->Index, DocumentBegin[NumDocuments - 1]->Index);
  CHECK(NumDocuments <= LBData.size());

  bool isInitialized = true;
  for (uint32_t I = 0; I < NumDocuments; I++) {
    const Document* Doc = DocumentBegin[I];
    if (CachedDocs.find(Doc->Index) == CachedDocs.end()) {
      isInitialized = false;
      break;
    }
  }

  if (isInitialized) {
    LOG_IF(verbose >= 2, "  already cached");
    return;
  }
  
  LOG_IF(verbose >= 2, "  init the cache for %d docs", NumDocuments);

  CachedDocs.clear();
  for (uint32_t I = 0; I < NumDocuments; I++) {
    Document* Doc = DocumentBegin[I];
    Doc->LocalIndex = I;
    CachedDocs.insert(Doc->Index);
  }

  for (uint32_t L = 0; L < NumDocuments; L++) {
    LBData[L][L] = 0;
    for (uint32_t R = L + 1; R < NumDocuments; R++) {
      const Document* DocL = DocumentBegin[L];
      const Document* DocR = DocumentBegin[R];
      LBData[DocL->LocalIndex][DocR->LocalIndex] = computeCrossings(DocL, DocR);
      LBData[DocR->LocalIndex][DocL->LocalIndex] = computeCrossings(DocR, DocL);
    }
  }
  LOG_IF(verbose >= 2, "  done init the cache for %d docs", NumDocuments);
}

template class BalancedPartitioning<uint8_t>;
template class BalancedPartitioning<uint16_t>;
template class BalancedPartitioning<uint32_t>;
template class BalancedPartitioning<uint64_t>;
