"""Thin wrapper around sentence-transformers for embedding text into vectors.

Uses the all-MiniLM-L6-v2 model — a small (22M param, 384-dim) model that
runs on CPU with no GPU or API key required.

Network behavior and privacy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No user code, queries, or file contents are EVER sent to HuggingFace or any
external service. All embedding computation happens locally. The only network
activity is the one-time model download on first use:

  First run (model not yet cached, ~80MB download):
    - Downloads model weight files from huggingface.co.
    - The only metadata sent in requests is a user-agent header with library
      versions (e.g. "hf_hub/1.6.0; python/3.12.3; torch/2.10.0") and a
      random session ID. No user-identifiable information is sent.
    - HF_HUB_DISABLE_TELEMETRY=1 is set, which disables any optional
      analytics the library might otherwise send.

  Subsequent runs (model already cached):
    - HF_HUB_OFFLINE=1 is set automatically (see code below), which prevents
      ALL network calls. The model loads entirely from the local cache at
      ~/.cache/huggingface/hub/.

This module is imported by both the MCP server (for search queries) and the
indexer (for bulk embedding of code chunks). The bulk embed() call shows a
progress bar; the single-query embed_query() does not.

The env vars below must be set before importing transformers or
sentence_transformers. They suppress warnings, progress bars, and telemetry.
Stray stderr output would interfere with the MCP server's JSON-RPC transport.
"""

import os
from pathlib import Path

os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_VERBOSITY", "error")
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

# If the model is already downloaded, force fully offline mode so that no
# network calls are made — not even the etag freshness check that
# huggingface_hub does by default. The first-ever run will still download
# normally because the cache dir won't exist yet. This MUST be set before
# importing sentence_transformers, which reads the env var at import time.
_MODEL_NAME = "all-MiniLM-L6-v2"
_HF_CACHE = Path(os.environ.get("HF_HOME", Path.home() / ".cache" / "huggingface")) / "hub"
if (_HF_CACHE / f"models--sentence-transformers--{_MODEL_NAME}").exists():
    os.environ.setdefault("HF_HUB_OFFLINE", "1")

import transformers  # noqa: E402 — must come after env vars are set
transformers.logging.disable_progress_bar()


class EmbeddingProvider:
    """Sentence-transformers embedder using all-MiniLM-L6-v2."""

    def __init__(self, model_name: str = _MODEL_NAME):
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
