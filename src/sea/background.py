"""Exact batched fingerprint comparisons and reusable Monte Carlo distributions.

Sampling order and floating-point addition order match statistics.significance's
scalar reference path. Only the calculation of fingerprint intersections changes.
"""
import os
from pathlib import Path
import secrets

import numpy as np
from rdkit import DataStructs


class FingerprintIndex:
    def __init__(self, pool, check=lambda: None):
        width = pool[0].fp.GetNumBits()
        if width % 64:
            raise ValueError("Batched fingerprints must have a multiple of 64 bits.")
        self.bits = np.empty((len(pool), width // 64), dtype=np.uint64)
        for i, molecule in enumerate(pool):
            if i % 4096 == 0:
                check()
            self.bits[i] = np.frombuffer(DataStructs.BitVectToBinaryText(molecule.fp), dtype=np.uint64)
        self.counts = np.bitwise_count(self.bits).sum(axis=1, dtype=np.int32)

    def scores(self, left, right, threshold, check):
        """Return one raw score per trial, preserving left-major pair order."""
        a, b = np.broadcast_arrays(left[:, :, None], right[:, None, :])
        ac, bc = self.counts[a], self.counts[b]
        maximum = np.maximum(ac, bc)
        upper = np.divide(np.minimum(ac, bc), maximum,
                          out=np.zeros(maximum.shape, dtype=float), where=maximum != 0)
        # min(on_bits) / max(on_bits) is an upper bound on Tanimoto.
        # Equality is retained. RDKit BulkTanimotoSimilarity assigns zero
        # (rather than one) to a pair of empty bit fingerprints.
        candidates = np.flatnonzero(upper.ravel() >= threshold)
        values = np.zeros(a.shape, dtype=float)
        for start in range(0, len(candidates), 8192):
            check()
            positions = candidates[start:start + 8192]
            ai, bi = a.flat[positions], b.flat[positions]
            common = np.bitwise_count(self.bits[ai] & self.bits[bi]).sum(axis=1, dtype=np.int32)
            union = self.counts[ai] + self.counts[bi] - common
            similarities = np.divide(common, union, out=np.zeros(len(union), dtype=float), where=union != 0)
            values.flat[positions] = np.where(similarities >= threshold, similarities, 0.)
        # Unlike a reduction/tree sum, accumulate uses the original sequential
        # double additions. Zero terms do not change the nonnegative score.
        return np.add.accumulate(values.reshape(len(left), -1), axis=1)[:, -1]


class NullCache(dict):
    """Worker memory cache backed by atomic, dataset-specific NumPy files."""
    def __init__(self, root, dataset_id):
        super().__init__()
        self.directory = Path(root) / "cache" / "null-v1" / dataset_id
        self.index = None

    def sample(self, pool, model, key, check):
        trials = model["empirical_trials"]
        path = self.directory / f"{key[0]}-{key[1]}.npy"
        check()
        try:
            values = np.load(path, allow_pickle=False)
            if (values.dtype == np.dtype("float64") and values.shape == (trials,)
                    and np.isfinite(values).all() and (values >= 0).all()
                    and (np.diff(values) >= 0).all()):
                return values
        except (OSError, ValueError, EOFError):
            pass
        if self.index is None:
            self.index = FingerprintIndex(pool, check)
        rng = np.random.default_rng(np.random.SeedSequence([model["seed"], *key]))
        scores = []
        block_size = max(1, min(64, 65536 // (key[0] * key[1])))
        for start in range(0, trials, block_size):
            check()
            count = min(block_size, trials - start)
            left = np.empty((count, key[0]), dtype=np.int64)
            right = np.empty((count, key[1]), dtype=np.int64)
            for i in range(count):
                # Keep the same alternating RNG calls as the scalar sampler.
                left[i] = rng.choice(len(pool), key[0], replace=False)
                right[i] = rng.choice(len(pool), key[1], replace=False)
            scores.extend(self.index.scores(left, right, model["threshold"], check))
        values = np.sort(scores)
        check()
        temporary = path.with_suffix(f".{os.getpid()}.{secrets.token_hex(4)}.tmp")
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            with temporary.open("wb") as stream:
                np.save(stream, values, allow_pickle=False)
            temporary.replace(path)
        except OSError:
            # A read-only/full cache directory must not invalidate an analysis.
            pass
        finally:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
        return values
