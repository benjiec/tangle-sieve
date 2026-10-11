"""Fold discovery across filesystem directories and nested ZIP archives."""
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


def ignored(name):
    return any(p == '__MACOSX' or p.startswith('._') for p in str(name).split('/'))


def model_number(name):
    if ignored(name):
        return None
    match = re.search(r'_model_(\d+)\.cif$', str(name))
    return int(match.group(1)) if match else None


@dataclass
class Fold:
    source_path: str
    fold_path: str
    models: list
    read: object


def discover(path, error):
    """Yield folds with readers valid until the generator advances; report local errors."""
    path = Path(path).resolve()

    def archive_folds(archive, source, prefix=''):
        groups = {}
        nested = []
        candidates = set()
        for info in archive.infolist():
            if info.is_dir() or ignored(info.filename):
                continue
            parent = PurePosixPath(info.filename).parent
            if (info.filename.endswith('_job_request.json') or
                    (info.filename.endswith('.cif') and 'templates' not in parent.parts)):
                candidates.add(str(parent))
            if model_number(info.filename) is not None:
                groups.setdefault(str(PurePosixPath(info.filename).parent), []).append(info)
            elif info.filename.lower().endswith('.zip'):
                nested.append(info)
        for directory in sorted(candidates - groups.keys()):
            error(f'{source}::{prefix}{directory}: fold candidate has no recognized model CIF files')
        for directory, infos in sorted(groups.items()):
            names = sorted(i.filename for i in infos)
            yield Fold(source, prefix + directory, names,
                       lambda name: archive.read(name).decode('utf-8'))
        for info in sorted(nested, key=lambda i: i.filename):
            try:
                with tempfile.SpooledTemporaryFile(max_size=16 * 1024 * 1024) as stream:
                    with archive.open(info) as member:
                        shutil.copyfileobj(member, stream)
                    stream.seek(0)
                    with zipfile.ZipFile(stream) as inner:
                        yield from archive_folds(inner, source, prefix + info.filename + '!/')
            except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as exc:
                error(f'{source}::{prefix}{info.filename}: {exc}')

    def walk(current):
        try:
            if current.is_dir():
                children = sorted(p for p in current.iterdir() if not ignored(p.name) and not p.is_symlink())
                models = [p for p in children if p.is_file() and model_number(p.name) is not None]
                if models:
                    yield Fold(str(current), '.', [p.name for p in models],
                               lambda name: (current / name).read_text(encoding='utf-8'))
                if not models and any(p.is_file() and (p.name.endswith('_job_request.json') or
                        (p.suffix == '.cif' and 'templates' not in current.parts)) for p in children):
                    error(f'{current}: fold candidate has no recognized model CIF files')
                for child in children:
                    if child.is_dir() or child.suffix.lower() == '.zip':
                        yield from walk(child)
            elif current.suffix.lower() == '.zip':
                with zipfile.ZipFile(current) as archive:
                    yield from archive_folds(archive, str(current))
            else:
                error(f'{current}: expected a directory or ZIP archive')
        except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as exc:
            error(f'{current}: {exc}')

    yield from walk(path)
