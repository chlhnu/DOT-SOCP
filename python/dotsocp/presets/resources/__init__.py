from pathlib import Path


def _resolve_resources_path(name: str) -> Path:
    """Resolve an resources path with the following priority:

    1. If `name` is an absolute or existing relative path: use it directly.
    2. Otherwise, look under this package's `resources/` directory:
         <this_file_dir>/<name>
    """
    # 1) Explicit / already valid relative path
    candidate = Path(name)
    if candidate.is_absolute() or candidate.exists():
        if candidate.exists():
            return candidate

    # 2) Local resources directory
    module_dir = Path(__file__).resolve().parent
    local_resource = module_dir / name
    if local_resource.exists():
        return local_resource

    # Not found anywhere
    search_paths = [
        str(candidate),
        str(local_resource),
    ]
    raise FileNotFoundError(
        f"Resource '{name}' not found. Tried:\n  - " + "\n  - ".join(search_paths)
    )
