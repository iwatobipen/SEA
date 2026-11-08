"""Resumable, atomic sampling checkpoints shared by fixed-budget profiles.

Unsorted samples and RNG state are committed together. A cancelled block is
discarded, so resumption is exactly equivalent to uninterrupted sampling.
"""
import copy
import fcntl
import hashlib
import json
from pathlib import Path
import secrets
import time
from zipfile import BadZipFile

import numpy as np

from .background import FingerprintIndex
from .precision import NULL_DEFINITION, SAMPLER_VERSION


class PrecisionCache(dict):
    def __init__(self, root, dataset_id, index_owner=None):
        super().__init__()
        self.directory = Path(root) / "cache" / "null-v2" / dataset_id
        self.index_owner = index_owner
        self.index = None
        self.on_progress = lambda *_: None

    def sample(self, pool, model, key, check):
        budget = model["empirical_trials"]
        identity = {"null": NULL_DEFINITION, "implementation": SAMPLER_VERSION,
                    "seed": model["seed"], "threshold": model["threshold"],
                    "sizes": list(key), "pool_size": len(pool)}
        token = json.dumps(identity, sort_keys=True)
        path = self.directory / (hashlib.sha256(token.encode()).hexdigest()[:20] + ".npz")
        # Atomic checkpoints are safe to read while another worker extends the
        # stream. Screening must not wait for a long refinement holding a lock.
        existing, _ = self._checkpoint(path, token)
        if len(existing) >= budget:
            check()
            self.on_progress(budget, budget)
            return np.sort(existing[:budget])
        lock = None
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            lock = path.with_suffix(".lock").open("a+b")
            while True:
                check()
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    time.sleep(.05)
        except OSError:
            if lock:
                lock.close()
            lock = None
        try:
            return self._sample(pool, model, key, check, path, token, budget, lock is not None)
        finally:
            if lock:
                lock.close()

    @staticmethod
    def _checkpoint(path, token):
        try:
            with np.load(path, allow_pickle=False) as archive:
                stored = archive["values"]
                state_text = str(archive["state"].item())
                digest = hashlib.sha256(stored.tobytes() + state_text.encode() + token.encode()).hexdigest()
                if (str(archive["identity"].item()) == token and str(archive["checksum"].item()) == digest
                        and stored.dtype == np.dtype("float64") and stored.ndim == 1
                        and np.isfinite(stored).all() and (stored >= 0).all()):
                    state = json.loads(state_text)
                    # Validate the RNG state before accepting any stored prefix.
                    np.random.default_rng().bit_generator.state = state
                    return stored.copy(), state
        except (OSError, ValueError, EOFError, KeyError, TypeError, BadZipFile):
            pass
        return np.empty(0, dtype=np.float64), None

    def _sample(self, pool, model, key, check, path, token, budget, writable):
        rng = np.random.default_rng(np.random.SeedSequence([model["seed"], *key]))
        existing, state = self._checkpoint(path, token)
        if state is not None:
            rng.bit_generator.state = state
        check()
        if len(existing) >= budget:
            self.on_progress(budget, budget)
            return np.sort(existing[:budget])
        values = np.empty(budget, dtype=np.float64)
        values[:len(existing)] = existing
        count = len(existing)
        stable_state = copy.deepcopy(rng.bit_generator.state)
        last_save = time.monotonic()

        def save():
            if not writable:
                return
            temporary = path.with_suffix(f".{secrets.token_hex(8)}.tmp")
            state_text = json.dumps(stable_state, sort_keys=True)
            checksum = hashlib.sha256(values[:count].tobytes() + state_text.encode() + token.encode()).hexdigest()
            try:
                with temporary.open("wb") as stream:
                    np.savez(stream, values=values[:count], state=state_text, identity=token, checksum=checksum)
                temporary.replace(path)
            except OSError:
                pass
            finally:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError:
                    pass

        try:
            if self.index is None:
                self.index = getattr(self.index_owner, "index", None)
                if self.index is None:
                    self.index = FingerprintIndex(pool, check)
                    if self.index_owner is not None:
                        self.index_owner.index = self.index
            block = max(1, min(64, 65536 // (key[0] * key[1])))
            while count < budget:
                check()
                size = min(block, budget - count)
                left = np.empty((size, key[0]), dtype=np.int64)
                right = np.empty((size, key[1]), dtype=np.int64)
                for i in range(size):
                    left[i] = rng.choice(len(pool), key[0], replace=False)
                    right[i] = rng.choice(len(pool), key[1], replace=False)
                batch = self.index.scores(left, right, model["threshold"], check)
                values[count:count + size] = batch
                count += size
                stable_state = copy.deepcopy(rng.bit_generator.state)
                if time.monotonic() - last_save >= 1.:
                    save()
                    self.on_progress(count, budget)
                    last_save = time.monotonic()
            self.on_progress(count, budget)
            return np.sort(values)
        finally:
            save()

    def distribution(self, pool, model, key, check):
        # In-memory arrays are keyed by the planned budget; different profiles
        # share the underlying unsorted checkpoint without sharing an incorrect
        # sorted prefix. The checkpoint preserves the original sample order.
        cache_key = (*key, model["empirical_trials"])
        if cache_key not in self:
            values = self.sample(pool, model, key, check)
            # Bound long-lived worker memory independently of library size.
            while self and sum(array.nbytes for array in self.values()) + values.nbytes > 64 * 1024**2:
                self.pop(next(iter(self)))
            self[cache_key] = values
        values = self.pop(cache_key)
        self[cache_key] = values
        return values
