"""Embedding provider using sentence-transformers (all-MiniLM-L6-v2).

Requires: pip install sentence-transformers
"""

import os

# Suppress noisy HuggingFace warnings about authentication
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")


class EmbeddingProvider:
    """Sentence-transformers embedder using all-MiniLM-L6-v2."""

    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        from sentence_transformers import SentenceTransformer
        self.model = SentenceTransformer(model_name)
        self.dim = self.model.get_sentence_embedding_dimension()

    def embed(self, texts: list[str]) -> list[list[float]]:
        """Embed a list of texts into vectors."""
        embeddings = self.model.encode(texts, show_progress_bar=True,
                                       batch_size=64)
        return embeddings.tolist()

    def embed_query(self, text: str) -> list[float]:
        """Embed a single query text."""
        return self.model.encode([text])[0].tolist()
