"""Cache global scaffold groups and a Morgan comparator for parameter validation.

Standardized structures are read from an immutable dataset. Chunk checkpoints
allow interrupted preparation to resume without touching application datasets.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from sea.data import load_dataset


def prepare_chunk(start, smiles, directory):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
    groups = np.empty(len(smiles), dtype='S64')
    fingerprints = np.empty((len(smiles), 256), dtype=np.uint8)
    for i, value in enumerate(smiles):
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f'Invalid standardized structure at row {start+i}')
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        groups[i] = hashlib.sha256(scaffold.encode()).hexdigest()
        fingerprints[i] = np.frombuffer(DataStructs.BitVectToBinaryText(generator.GetFingerprint(mol)), dtype=np.uint8)
    path = Path(directory) / f'{start}.npz'
    temporary = path.with_suffix('.tmp')
    with temporary.open('wb') as stream:
        np.savez(stream, groups=groups, fingerprints=fingerprints)
    temporary.replace(path)
    return start, len(smiles)


def prepare(dataset_id, directory, workers):
    started = time.perf_counter()
    data, structures = load_dataset(dataset_id)
    smiles = list(structures)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    metadata = {'dataset_id': dataset_id, 'molecules': len(smiles), 'rdkit_version': rdBase.rdkitVersion,
                'smiles_order_sha256': hashlib.sha256('\n'.join(smiles).encode()).hexdigest(),
                'scaffold': 'MurckoScaffoldSmiles(includeChirality=False); SHA256; global grouping',
                'comparator': {'kind': 'Morgan', 'radius': 2, 'bits': 2048, 'include_chirality': False}}
    identity = directory / 'identity.json'
    if identity.exists() and json.loads(identity.read_text()) != metadata:
        raise ValueError('Cache identity differs. Choose a new directory.')
    identity.write_text(json.dumps(metadata, indent=2) + '\n')
    if (directory / 'complete.json').exists():
        print('Preparation cache already complete.', flush=True)
        return
    chunk_size = 4096
    starts = list(range(0, len(smiles), chunk_size))
    missing = [start for start in starts if not (directory / f'{start}.npz').exists()]
    print(f'Preparing {len(smiles):,} structures; {len(missing)} chunks remain.', flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(prepare_chunk, start, smiles[start:start+chunk_size], str(directory)) for start in missing]
        for i, future in enumerate(as_completed(futures), 1):
            future.result()
            if i % 5 == 0 or i == len(futures):
                print(f'Chunks complete: {len(starts)-len(missing)+i}/{len(starts)}; {time.perf_counter()-started:.1f} s', flush=True)
    groups = np.lib.format.open_memmap(directory / 'groups.tmp.npy', mode='w+', dtype='S64', shape=(len(smiles),))
    fps = np.lib.format.open_memmap(directory / 'morgan.tmp.npy', mode='w+', dtype=np.uint8, shape=(len(smiles), 256))
    for start in starts:
        with np.load(directory / f'{start}.npz', allow_pickle=False) as archive:
            count = len(archive['groups'])
            groups[start:start+count] = archive['groups']
            fps[start:start+count] = archive['fingerprints']
    groups.flush()
    fps.flush()
    del groups, fps
    (directory / 'groups.tmp.npy').replace(directory / 'groups.npy')
    (directory / 'morgan.tmp.npy').replace(directory / 'morgan.npy')
    (directory / 'complete.json').write_text(json.dumps({**metadata, 'elapsed_seconds': time.perf_counter()-started}, indent=2) + '\n')
    print(f'Preparation complete: {time.perf_counter()-started:.1f} s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-id', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('Use 1 to 8 preparation workers.')
    prepare(args.dataset_id, args.cache, args.workers)
