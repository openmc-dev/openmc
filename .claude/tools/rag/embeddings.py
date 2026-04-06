"""Thin wrapper around sentence-transformers for embedding text into vectors.

Uses the all-MiniLM-L6-v2 model — a small (22M param, 384-dim) model that
runs on CPU with no GPU or API key required.

Network behavior and privacy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No user code, queries, or file contents are EVER sent to HuggingFace or any
external service. All embedding computation happens locally. The only network
activity is the one-time model download on first use:

  First run (model not yet cached, ~80MB download):
    - Downloads model weight files from huggingface.co. This is a standard
      HTTP file download, similar to pip installing a package.
    - The only metadata sent in these requests is an HTTP user-agent header
      containing library version numbers (e.g. "hf_hub/1.6.0;
      python/3.12.3; torch/2.10.0"). No filenames, file contents, queries,
      or any user-identifiable information is sent.
    - The huggingface_hub library has an optional feature where it can report
      anonymous library usage statistics (just version numbers, not user
      data) back to HuggingFace. We disable this by setting
      HF_HUB_DISABLE_TELEMETRY=1.

  Subsequent runs (model already cached):
    - We set HF_HUB_OFFLINE=1 automatically (see _set_offline_if_cached()
      below), which prevents ALL network calls. The model loads entirely
      from the local cache at ~/.cache/huggingface/hub/. Zero bytes leave
      the machine.

How the model is downloaded
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The SentenceTransformer() constructor (called in __init__ below) handles
the download automatically on first use. It calls into the huggingface_hub
library, which downloads the model files from:

    https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2

The files are saved to ~/.cache/huggingface/hub/ and reused on subsequent
runs. We pass token=False to ensure no authentication token is sent.

This module is imported by both the MCP server (for search queries) and the
indexer (for bulk embedding of code chunks). The bulk embed() call shows a
progress bar; the single-query embed_query() does not.

The env vars below must be set before importing transformers or
sentence_transformers. They suppress warnings and progress bars that these
libraries emit by default. Stray stderr output would interfere with the MCP
server's JSON-RPC transport.
"""

import os
from pathlib import Path

MODEL_NAME = "all-MiniLM-L6-v2"

# These env vars control logging behavior in the HuggingFace libraries.
# They must be set before the libraries are imported.
os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")   # suppress warnings
os.environ.setdefault("HF_HUB_VERBOSITY", "error")        # suppress warnings
os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")   # suppress threading warning
# Disable anonymous library usage statistics (version numbers only, not user
# data — but we disable it anyway as a matter of policy).
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")


def _set_offline_if_cached():
    """If the model has already been downloaded, tell huggingface_hub to
    skip all network calls by setting HF_HUB_OFFLINE=1.

    Without this, huggingface_hub makes an HTTP request to huggingface.co
    on every load to check if the cached model is still up to date — even
    though the model never changes. Setting HF_HUB_OFFLINE=1 prevents this.

    This must run before sentence_transformers is imported, because the
    library reads the env var at import time.
    """
    # HuggingFace caches downloaded models under ~/.cache/huggingface/hub/
    # in directories named like "models--sentence-transformers--all-MiniLM-L6-v2".
    # The HF_HOME env var can override the base cache location.
    hf_home = os.environ.get("HF_HOME")
    if hf_home:
        cache_dir = Path(hf_home) / "hub"
    else:
        cache_dir = Path.home() / ".cache" / "huggingface" / "hub"

    model_dir = cache_dir / f"models--sentence-transformers--{MODEL_NAME}"
    if model_dir.exists():
        os.environ.setdefault("HF_HUB_OFFLINE", "1")


_set_offline_if_cached()

# This import must come after the env vars above are set, because the
# transformers library reads them at import time.
import transformers
transformers.logging.disable_progress_bar()


class EmbeddingProvider:
    """Sentence-transformers embedder using all-MiniLM-L6-v2."""

    def __init__(self, model_name: str = MODEL_NAME):
        from sentence_transformers import SentenceTransformer

        # This constructor loads the model from the local cache. If the model
        # has not been downloaded yet, it downloads it from huggingface.co
        # (~80MB, one-time). token=False ensures no auth token is sent.
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
