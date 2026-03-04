"""Embedding provider with auto-detection fallback chain.

1. sentence-transformers (all-MiniLM-L6-v2) - good quality, ~80MB model
2. TF-IDF + SVD - zero downloads, decent for code identifiers
"""

import os
import sys
from abc import ABC, abstractmethod

# Suppress noisy HuggingFace warnings about authentication
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")


class EmbeddingProvider(ABC):
    """Abstract base for embedding providers."""

    dim: int = 0

    @abstractmethod
    def embed(self, texts: list[str]) -> list[list[float]]:
        """Embed a list of texts into vectors."""
        ...

    @abstractmethod
    def embed_query(self, text: str) -> list[float]:
        """Embed a single query text."""
        ...

    @staticmethod
    def create(corpus_texts: list[str] | None = None) -> "EmbeddingProvider":
        """Auto-detect best available embedding backend.

        Args:
            corpus_texts: For TF-IDF fallback, the full corpus to fit on.
                          Not needed for sentence-transformers.
        """
        # Try sentence-transformers first
        try:
            return SentenceTransformerProvider()
        except (ImportError, Exception) as e:
            print(f"  sentence-transformers unavailable: {e}", file=sys.stderr)

        # Fall back to TF-IDF
        if corpus_texts is None:
            raise RuntimeError(
                "No embedding provider available. Install sentence-transformers "
                "or provide corpus_texts for TF-IDF fallback."
            )
        print("  Using TF-IDF fallback embeddings", file=sys.stderr)
        return TfidfProvider(corpus_texts)


class SentenceTransformerProvider(EmbeddingProvider):
    """sentence-transformers with all-MiniLM-L6-v2."""

    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        from sentence_transformers import SentenceTransformer
        self.model = SentenceTransformer(model_name)
        self.dim = self.model.get_sentence_embedding_dimension()

    def embed(self, texts: list[str]) -> list[list[float]]:
        embeddings = self.model.encode(texts, show_progress_bar=True,
                                       batch_size=64)
        return embeddings.tolist()

    def embed_query(self, text: str) -> list[float]:
        return self.model.encode([text])[0].tolist()


class TfidfProvider(EmbeddingProvider):
    """TF-IDF vectors projected to dense via SVD. No model download needed."""

    def __init__(self, corpus_texts: list[str], dim: int = 256):
        from sklearn.decomposition import TruncatedSVD
        from sklearn.feature_extraction.text import TfidfVectorizer

        self.dim = dim
        self.vectorizer = TfidfVectorizer(
            max_features=10000,
            sublinear_tf=True,
            token_pattern=r"(?u)\b[a-zA-Z_][a-zA-Z0-9_]{2,}\b",  # Code identifiers
        )
        tfidf_matrix = self.vectorizer.fit_transform(corpus_texts)

        # Project to dense using SVD
        actual_dim = min(dim, tfidf_matrix.shape[1] - 1, tfidf_matrix.shape[0] - 1)
        self.svd = TruncatedSVD(n_components=actual_dim)
        self.svd.fit(tfidf_matrix)
        self.dim = actual_dim

    def embed(self, texts: list[str]) -> list[list[float]]:
        tfidf = self.vectorizer.transform(texts)
        dense = self.svd.transform(tfidf)
        return dense.tolist()

    def embed_query(self, text: str) -> list[float]:
        tfidf = self.vectorizer.transform([text])
        dense = self.svd.transform(tfidf)
        return dense[0].tolist()
