"""Thin wrapper around sentence-transformers for embedding text into vectors.

Uses the all-MiniLM-L6-v2 model — a small (22M param, 384-dim) model that
runs on CPU with no GPU or API key required. The model weights are downloaded
once from Hugging Face on first use and cached locally (~80MB).

This module is imported by both the MCP server (for search queries) and the
indexer (for bulk embedding of code chunks). The bulk embed() call shows a
progress bar; the single-query embed_query() does not.

The env vars and logging suppression below are necessary because HuggingFace
libraries emit warnings, progress bars, and telemetry pings by default. When
running under the MCP server, any stray output can interfere with the JSON-RPC
transport. These must be set before importing transformers/sentence_transformers.
"""

import os

os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

import transformers  # noqa: E402 — must come after env vars are set
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
