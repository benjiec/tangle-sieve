"""Fold discovery across filesystem directories and nested ZIP archives."""
import re
import json
import shutil
import struct
import tempfile
import zipfile
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
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
    member_names: tuple = ()
    timestamp: str = ''


def request_name(model_name):
    return re.sub(r'model_\d+\.cif$', 'job_request.json', str(model_name))


def read_seed(read, models, names=None):
    requests = {request_name(name) for name in models}
    if len(requests) != 1:
        raise ValueError('fold must have exactly one matching job-request file')
    name = next(iter(requests))
    if names is not None and list(names).count(name) != 1:
        raise ValueError(f'expected exactly one seed job-request file: {name}')
    try:
        request = json.loads(read(name))
    except (OSError, KeyError) as exc:
        raise ValueError(f'cannot read seed job-request file: {name}') from exc
    if isinstance(request, list):
        if len(request) != 1:
            raise ValueError(f'expected exactly one job in {name}')
        request = request[0]
    if not isinstance(request, dict):
        raise ValueError(f'invalid job request: {name}')
    seeds = request.get('modelSeeds')
    if not isinstance(seeds, list) or len(seeds) != 1:
        raise ValueError(f'expected exactly one modelSeeds value in {name}')
    seed = seeds[0]
    if isinstance(seed, bool) or not (isinstance(seed, int) or
                                    isinstance(seed, str) and re.fullmatch(r'[0-9]+', seed)):
        raise ValueError(f'invalid seed in {name}: {seed!r}')
    seed = int(seed)
    if seed < 0:
        raise ValueError(f'invalid negative seed in {name}')
    return seed


def zip_timestamp(infos):
    """Latest model mtime: prefer UTC extended timestamps, else unzoned ZIP times."""
    utc_times = []
    for info in infos:
        extra = info.extra
        offset = 0
        epoch = None
        while offset + 4 <= len(extra):
            tag, size = struct.unpack_from('<HH', extra, offset)
            offset += 4
            payload = extra[offset:offset + size]
            offset += size
            if tag == 0x5455 and len(payload) >= 5 and payload[0] & 1:
                epoch = struct.unpack_from('<I', payload, 1)[0]
                break
        if epoch is not None:
            utc_times.append(datetime.fromtimestamp(epoch, timezone.utc))
    if len(utc_times) == len(infos):
        return max(utc_times).isoformat(timespec='seconds').replace('+00:00', 'Z')
    # DOS ZIP timestamps do not encode a timezone; never invent one.
    return max(datetime(*info.date_time) for info in infos).isoformat(timespec='seconds')


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
                       lambda name: archive.read(name).decode('utf-8'),
                       member_names=tuple(archive.namelist()), timestamp=zip_timestamp(infos))
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
                               lambda name: (current / name).read_text(encoding='utf-8'),
                               member_names=tuple(p.name for p in children if p.is_file()),
                               timestamp=datetime.fromtimestamp(max(p.stat().st_mtime for p in models),
                                                                timezone.utc).isoformat(timespec='seconds').replace('+00:00', 'Z'))
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


def fold_locator(source_path, fold_path):
    """Serialize the filesystem source and optional archive member location."""
    if fold_path == '.' and Path(source_path).suffix.lower() != '.zip':
        return source_path
    return f'{source_path}::{fold_path}'




@contextmanager
def open_fold(locator):
    """Open exactly one archive fold, without visiting unrelated nested archives."""
    source, separator, member_path = str(locator).partition('::')
    if not separator or not source or not member_path:
        raise ValueError('fold locator requires SOURCE.zip::FOLD_PATH')
    parts = member_path.split('!/')
    if any(not part for part in parts):
        raise ValueError('fold locator has an empty archive or directory component')
    with ExitStack() as stack:
        archive = stack.enter_context(zipfile.ZipFile(source))

        def read_unique(name):
            if archive.namelist().count(name) != 1:
                raise ValueError(f'expected exactly one archive member: {name}')
            return archive.read(name)

        for member in parts[:-1]:
            if archive.namelist().count(member) != 1:
                raise ValueError(f'expected exactly one nested archive: {member}')
            stream = stack.enter_context(tempfile.SpooledTemporaryFile(max_size=16 * 1024 * 1024))
            with archive.open(member) as inner:
                shutil.copyfileobj(inner, stream)
            stream.seek(0)
            archive = stack.enter_context(zipfile.ZipFile(stream))
        directory = parts[-1]
        names = [name for name in archive.namelist()
                 if str(PurePosixPath(name).parent) == directory and model_number(name) is not None]
        if not names:
            raise ValueError(f'no model CIF files at fold locator: {locator}')
        yield Fold(str(Path(source).resolve()), member_path, sorted(names),
                   lambda name: read_unique(name).decode('utf-8'), tuple(archive.namelist()))
