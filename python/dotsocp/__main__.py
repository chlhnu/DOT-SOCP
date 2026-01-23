from __future__ import annotations

try:
    from .cli import main
except ImportError:
    from dotsocp.cli import main


if __name__ == "__main__":
    main()
