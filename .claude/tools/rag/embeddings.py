"""Embedding provider using sentence-transformers (all-MiniLM-L6-v2).

Requires: pip install sentence-transformers
"""

import os

# Official HuggingFace/transformers knobs for quiet operation.
# TRANSFORMERS_VERBOSITY: controls transformers' own logging (load reports, etc.)
# HF_HUB_DISABLE_TELEMETRY: don't phone home
# TOKENIZERS_PARALLELISM: avoids fork-safety warning
# All must be set before importing transformers.
os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

import transformers
transformers.logging.disable_progress_bar()


class EmbeddingProvider:
    """Sentence-transformers embedder using all-MiniLM-L6-v2."""

    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        from sentence_transformers import SentenceTransformer
        self.model = SentenceTransformer(model_name, token=False)
        self.dim = self.model.get_sentence_embedding_dimension()

    def embed(self, texts: list[str]) -> list[list[float]]:
        """Embed a list of texts into vectors."""
        embeddings = self.model.encode(texts, show_progress_bar=True,
                                       batch_size=64)
        return embeddings.tolist()

    def embed_query(self, text: str) -> list[float]:
        """Embed a single query text."""
        return self.model.encode([text])[0].tolist()
