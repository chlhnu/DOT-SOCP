from __future__ import annotations

from typing import Iterable, List, Optional
import logging
from dotsocp.configs.icons import LOG_ICON
from dotsocp.solvers.algorithms.utils.admm_utils import RunHistoryManager


class MultilevelRunHistoryManager:
    """
    Lightweight wrapper that holds per-level RunHistoryManager objects and can
    print a combined summary. Intended as a minimal adapter around the existing
    per-level histories.
    """

    def __init__(self, all_time: float, histories: Iterable[RunHistoryManager] | None = None) -> None:
        self.all_time = all_time
        self.histories: List[RunHistoryManager] = list(histories or [])

    def append(self, history: RunHistoryManager) -> None:
        self.histories.append(history)

    def __iter__(self):
        return iter(self.histories)

    def __len__(self) -> int:
        return len(self.histories)

    def last(self) -> RunHistoryManager | None:
        return self.histories[-1] if self.histories else None
    
    def print(self) -> None:
        last_history = self.last()

        if last_history is None:
            return
        
        logging.info("%s  Last Level Summary", LOG_ICON.get("stop", "[history]"))
        last_history.print_end_history()
        last_history.print_steps_time()

        level_times = [round(float(h.running_time), 2) for h in self.histories]
        level_iters = [int(h.kkt_iteration[-1]) for h in self.histories]

        label_width = 19
        logging.info(
            f"{LOG_ICON.get('stop', '[history]')}  Multilevel Summary\n"
            f"{"Iteration of levels".ljust(label_width)}: {level_iters}\n"
            f"{"Time of levels".ljust(label_width)}: {level_times} sec\n"
            f"{"Multilevel all time".ljust(label_width)}: {self.all_time:.2f} sec\n"
            + '=' * 60 + '\n'
        )
    
    def plot(self, filename: Optional[str] = None) -> None:
        last_history = self.last()

        if last_history is None:
            logging.info(f"{LOG_ICON.get('stop', '[history]')}  No history to plot KKT error curves.")
            return
        
        logging.info(f"{LOG_ICON.get('stop', '[history]')}  Plotting KKT error curves of the last level...")
        last_history.plot_kkt_errors(filename=filename)
    