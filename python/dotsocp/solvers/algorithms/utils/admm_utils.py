import sys
import time

import numpy as np
import logging
from typing import List, Any, Union, Optional
from matplotlib import pyplot as plt
from contextlib import contextmanager

from numpy import ndarray
from math import log10

from tqdm import tqdm as tqdm_std # text-mode tqdm (works in terminals and notebooks)
from tqdm.auto import tqdm as tqdm_auto # auto-selecting tqdm

def _in_notebook() -> bool:
    try:
        from IPython import get_ipython  # type: ignore

        ip = get_ipython()
        if ip is None:
            return False

        return ip.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


# Only show tqdm progress bar when output to an interactive terminal
def _smart_tqdm(*args, **kwargs):
    # Respect explicit user choice if provided.
    if "disable" not in kwargs and (not sys.stdout.isatty()) and (not _in_notebook()):
        kwargs['disable'] = True

    # In notebooks, prefer the text-mode progress bar for predictable formatting.
    tqdm_impl = tqdm_std if _in_notebook() else tqdm_auto
    return tqdm_impl(*args, **kwargs)


class PenaltyParamManager:
    """A history-based penalty parameter manager for ADMM.
    """

    def __init__(self):
        self.last_it = -float('inf')
        self._sigma_upper_bound = 10 ** 3
        self._sigma_lower_bound = 10 ** (-3)

    # ==== Adjust penalty parameter
    def is_to_adjust(self, current_it):
        """Determine whether to adjust the penalty parameter based on current iteration and
        iterations passed since last adjustment.

        Parameters:
            current_it: Current iteration number

        Returns:
            bool: True if penalty parameter should be adjusted, False otherwise
        """
        passed_it = current_it - self.last_it

        # Adjustment frequency increases with iteration count
        if (current_it < 20 and passed_it >= 3) or \
                (current_it < 50 and passed_it >= 6) or \
                (current_it < 100 and passed_it >= 10) or \
                (current_it < 200 and passed_it >= 15) or \
                (current_it < 500 and passed_it >= 25) or \
                (passed_it >= 40):
            self.last_it = current_it
            return True

        return False

    def get_updated_value(self, sigma, prim_dual_ratio):
        """Adjust ADMM penalty parameter according to primal/dual ratio
        """

        # Update with safeguard
        adjust_factor = self.__get_adjust_factor(prim_dual_ratio)
        sigma = max(min(sigma * adjust_factor, self._sigma_upper_bound), self._sigma_lower_bound)

        return sigma

    @staticmethod
    def __get_adjust_factor(prim_dual_ratio):
        """Get adjusting factor of ADMM penalty parameter
        """

        # Normalize ratio to always be >= 1 for table lookup
        if prim_dual_ratio < 1.:
            lookup_ratio = 1. / prim_dual_ratio
            invert = True
        else:
            lookup_ratio = prim_dual_ratio
            invert = False

        # Find the appropriate factor
        match lookup_ratio:
            case xi if xi > 50:     factor = 2.00
            case xi if xi > 40:     factor = 1.80
            case xi if xi > 20:     factor = 1.60
            case xi if xi > 10:     factor = 1.40
            case xi if xi > 5.0:    factor = 1.35
            case xi if xi > 3.33:   factor = 1.32
            case xi if xi > 2.5:    factor = 1.28
            case xi if xi > 2.0:    factor = 1.26
            case xi if xi > 1.5:    factor = 1.20
            case xi if xi > 1.2:    factor = 1.15
            case xi if xi > 1.1:    factor = 1.10
            case _:                 factor = 1.00

        # Recover factor
        return 1. / factor if invert else factor


class RunHistoryManager:
    """A class to help record Optimization Algorithm running history and to help visualize them
    """

    def __init__(self, max_record_numbers: int, kkt_labels: List[str], name: str,
                       kkt_short_labels: Optional[List[str]] = None, use_linear_progress: bool = False):
        assert kkt_short_labels is None or len(kkt_short_labels) == len(kkt_labels), \
            "kkt_short_labels must be None or have the same length as kkt_labels."

        # Private value
        self._kkt_num = 0
        self._max_num = max_record_numbers
        self._unknown_value = np.inf
        self._start_time_global = self._unknown_value # global clock to time the algorithm and the kkt errors
        self._tol_progress = None # show progress of solution accuracy
        self._is_ended = True # Whether the end() function has been called after start()
        self._use_linear_progress = use_linear_progress # Whether to use linear progress bar
        self._separate_symbol = lambda s: f"{'----'} {s} ".ljust(42, '-')
        self._sub_separate_symbol = '-' * 42

        # Public value
        self.kkt_entry_num      = len(kkt_labels) # The number of kkt constraint
        self.kkt_errors         = np.ones((max_record_numbers, self.kkt_entry_num)) * self._unknown_value
        self.kkt_iteration      = np.ones(max_record_numbers) * self._unknown_value # The number of iteration which producing corresponding kkt errors
        self.kkt_time           = np.ones(max_record_numbers) * self._unknown_value # The time producing corresponding kkt errors
        self.kkt_labels         = kkt_labels
        self.kkt_short_labels   = kkt_short_labels if kkt_short_labels is not None else kkt_labels
        self.name               = name
        self.running_time       = self._unknown_value
        self.last_record_it     = -1 # The last iteration number that has been recorded

        self.steps_time: dict   = {} # Total running time of each step
        self.history: dict      = {} # Other running history besides kkt errors, which will be filled in self.record()

        logging.basicConfig(level=logging.INFO, format='%(message)s')

    # ================
    # Time for the whole algorithm
    def start(self):
        """Start running algorithm
        """

        self._start_time_global = time.perf_counter()
        self._is_ended = False

    def end(self):
        """Trim running history: truncate the redundant pre-allocated data
        """

        self.running_time = time.perf_counter() - self._start_time_global

        # Truncation
        self.kkt_errors = self.kkt_errors[:self._kkt_num, :]
        self.kkt_iteration = self.kkt_iteration[:self._kkt_num]
        self.kkt_time = self.kkt_time[:self._kkt_num]

        for key in self.history.keys():
            self.history[key] = self.history[key][:self._kkt_num]

        # Close
        if self._tol_progress is not None:
            self._tol_progress.close()
            # logging.info(self._separate_symbol("Finish performing"))
            sys.stdout.flush()
            time.sleep(0.1)
        
        self._is_ended = True
    
    def is_ended(self):
        """Check whether the end() function has been called after start()
        """
        return self._is_ended
    
    def get_running_time(self):
        """Gets the time that the algorithm has been run
        """
        return time.perf_counter() - self._start_time_global

    # ================
    @contextmanager
    def timer(self, tag):
        start_time = time.perf_counter()
        yield
        if tag in self.steps_time:
            self.steps_time[tag] += time.perf_counter() - start_time
        else:
            self.steps_time[tag] = time.perf_counter() - start_time

    def __transfer_tol_to_progress(self, tol):
        """Transfer tolerance (0 < tolerance < 1) to the integer type progress
        """
        if hasattr(self, '_use_linear_progress') and self._use_linear_progress:
            return self.__transfer_tol_to_progress_linear(tol)
        return self.__transfer_tol_to_progress_sublinear(tol)

    @staticmethod
    def __transfer_tol_to_progress_linear(tol):
        """Transfer tolerance (0 < tolerance < 1) to the integer type progress (linear rate)
        """
        return round(1000.0 * log10(1.0 / tol))
    
    @staticmethod
    def __transfer_tol_to_progress_sublinear(tol):
        """Transfer tolerance (0 < tolerance < 1) to the integer type progress (sublinear rate)
        """
        return round(1000.0 * (1.0 / tol) ** 0.5)

    def __format_condition_names(self, conditions):
        """Format condition indices/names for display
        
        Args:
            conditions: List of condition indices (int) or names (str)
            
        Returns:
            str: Formatted string of condition names
        """
        if not conditions:
            return "None"
        
        if isinstance(conditions[0], int):
            # Convert indices to short labels
            names = [self.kkt_short_labels[i] if i < len(self.kkt_short_labels) else f"KKT-{i}" 
                    for i in conditions]
        else:
            # Already names
            names = conditions
        
        # Limit display length
        if len(names) <= 2:
            return ", ".join(names)
        elif len(names) <= 4:
            return ", ".join(names[:2]) + f" + {len(names)-2} more"
        else:
            return f"{names[0]}, {names[1]} + {len(names)-2} others"

    def __get_condition_index(self, condition_name):
        """Get condition index from name"""
        try:
            return self.kkt_short_labels.index(condition_name)
        except ValueError:
            try:
                return self.kkt_labels.index(condition_name)
            except ValueError:
                return -1
            
    def __create_progress_bar(self):
        kwargs = {
            "total": self.__transfer_tol_to_progress(self._target_tol),
            "desc": f"Tol={self._target_tol:.2e}",
            "leave": True,
            "bar_format": "[{desc}{postfix}]|{bar}|{percentage:4.1f}%",
            # "dynamic_ncols": True,
            "ncols": 120,
        }

        return _smart_tqdm(**kwargs)

    def create_tol_progress(self, target_tol):
        """Create the progress object to show progress during iteration"""
        # logging.info(self._separate_symbol("Starting to perform ..."))
        
        # Store target tolerance for later use
        self._target_tol = target_tol
        self._tol_progress = self.__create_progress_bar()
        logging.debug("---- Iteration Start ".ljust(42, '-'))

    def show_tol_progress(self, current_it, current_err, active_idx=None, converged_idx=None):
        """Show the progress of solution accuracy w.r.t. target tolerance
        
        Args:
            current_it: Current iteration number
            current_err: Current maximum error value
            active_idx: List of indices or names of currently active KKT conditions
            converged_idx: List of indices or names of recently converged conditions
        """
        
        # Handle converged conditions - close current progress bar and start new one
        if converged_idx and len(converged_idx) > 0:
            if self._tol_progress is not None:
                # Show convergence message
                converged_names = self.__format_condition_names(converged_idx)
                self._tol_progress.set_postfix_str(f"Converged: {converged_names}")
                self._tol_progress.close()
                
                logging.info(f"Conditions converged at iteration {current_it}: {converged_names}\n")
            
            # Start new progress bar for remaining conditions
            converged_indices = converged_idx if (converged_idx and isinstance(converged_idx[0], int)) else \
                               [self.__get_condition_index(name) for name in converged_idx] if converged_idx else []
            
            # Get previously converged conditions (stored in instance)
            if not hasattr(self, '_converged_conditions'):
                self._converged_conditions = set()
            self._converged_conditions.update(converged_indices)
            
            remaining_conditions = [i for i in range(self.kkt_entry_num) 
                                  if i not in self._converged_conditions]
            
            if len(remaining_conditions) > 0:
                self._tol_progress = self.__create_progress_bar()
            else:
                logging.info("All KKT conditions have converged!")
                return
        
        # Update progress bar
        if self._tol_progress is not None:
            self._tol_progress.n = self.__transfer_tol_to_progress(max(current_err, self._target_tol))
            try:
                self._tol_progress.refresh()
            except Exception:
                pass

            # Calculate elapsed time since algorithm start
            elapsed_time = time.perf_counter() - self._start_time_global
            
            # Build postfix string manually to control order
            is_over_one_hour = (elapsed_time // 3600 >= 1)
            time_str = time.strftime('%M:%S', time.gmtime(elapsed_time)) if not is_over_one_hour else \
                time.strftime('%H:%M:%S', time.gmtime(elapsed_time))
            postfix_parts = [
                f"Acc: {current_err:.2e}",
                f"Time: {time_str}",
                f"Iter: {current_it} ({elapsed_time / (current_it + 1):.4f} sec/it)",
            ]

            if active_idx is not None:
                postfix_parts.append(f"Checking: {self.__format_condition_names(active_idx)}")

            self._tol_progress.set_postfix_str(", ".join(postfix_parts))
            self.__verbose_logging()

    def __verbose_logging(self):
        """Log the detailed information during iteration
        """
        idx = self._kkt_num - 1
        msg_it = f"{self.kkt_iteration[idx]:4.0f}"
        msg_kkt = " ".join([f"{error:6.2e}" for error in self.kkt_errors[idx, :]])
        logging.debug(f"Iteration: {msg_it} - KKT: {msg_kkt}")

    # ================
    # Record running history
    def record(self, current_it: Optional[int] = None,
                     kkt_errors: Optional[Union[List, ndarray]] = None,
                     history: Optional[dict[str, Any]] = None):
        """Record kkt errors and their corresponding iteration number and time
        """

        if kkt_errors is None or current_it is None:
            raise ValueError("Argument `kkt_errors` or `current_it` must be provided.")
        
        if current_it < self.last_record_it:
            raise ValueError(f"Current iteration {current_it} is smaller than last recorded iteration {self.last_record_it}. "
                             f"Please check the input arguments.")
        elif current_it == self.last_record_it:
            # If the current iteration is the same as the last recorded iteration,
            # we will overwrite the last record
            self._kkt_num -= 1
             

        if self._kkt_num < self._max_num:
            self.last_record_it = current_it

            self.kkt_errors[self._kkt_num, :] = kkt_errors
            self.kkt_iteration[self._kkt_num] = current_it
            self.kkt_time[self._kkt_num] = time.perf_counter() - self._start_time_global
        else:
            raise ValueError(f"Current recorded items number {self._kkt_num} is bigger than max valid number {self._max_num}."
                             f"There is no redundant space to store the running history.")

        if history is not None:
            for key, val in history.items():
                if key not in self.history:
                    self.history[key] = np.ones_like(self.kkt_iteration) * self._unknown_value
                self.history[key][self._kkt_num] = val

        self._kkt_num += 1
    
    def get_current_kkt_errors(self):
        """Get the current kkt errors
        """
        if self._kkt_num == 0:
            # raise ValueError("No kkt errors recorded yet.")
            return np.ones(self.kkt_errors.shape[1]) * self._unknown_value
        else:
            return self.kkt_errors[self._kkt_num - 1, :]

    # ================
    # Plot/Print running history
    def plot_kkt_errors(self, filename: Optional[str] = None, is_show_when_save: bool = False, x_axis: str = 'iteration',
                        title: Optional[str] = None, x_label: Optional[str] = None, y_label: Optional[str] = None):
        """Plot kkt error curves

        Usage
            1. plot_kkt_errors() to *show* kkt error curves
            2. plot_kkt_errors(file_name="name.ext") to *save* figure as a picture file (*not show* the figure)
            3. plot_kkt_errors(file_name="name.ext", is_show_when_save=True) to *show* and *save* the figure of kkt error curves
            4. plot_kkt_errors(x_axis=...)
                x_axis="iteration"(default), to set the x-axis label as iteration numbers
                x_axis="time", to set the x-axis label as iteration time
            5. fig = plot_kkt_errors(...) to get the figure handle
        """

        # Choose type of x-axis
        if x_axis == 'iteration':
            x_data = self.kkt_iteration
            x_default_label = "Iteration numbers"
        elif x_axis == 'time':
            x_data = self.kkt_time
            x_default_label = "Iteration time [seconds]"
        else:
            raise ValueError(f"x_axis {x_axis} is not supported.")

        # Plot
        fig = plt.figure()

        for n in range(self.kkt_entry_num):
            kkt_errors = self.kkt_errors[:, n]
            kkt_errors[kkt_errors < 10**(-10)] = 0.0 # If error < 1e-10, we consider it has been exactly solved
            plt.semilogy(x_data, kkt_errors, label=self.kkt_short_labels[n])

        plt.title(title if isinstance(title, str) else self.name)
        plt.xlabel(x_label if isinstance(x_label, str) else x_default_label)
        plt.ylabel(y_label if isinstance(y_label, str) else "Karush–Kuhn–Tucker errors")
        plt.legend()

        # Show / Save
        if isinstance(filename, str):
            if is_show_when_save:
                fig.show()

            format_fig = self._get_saved_fig_opts_matplotlib()
            fig.savefig(filename, **format_fig)
            plt.close(fig)
        else:
            fig.show()

        # plt.close(fig)

    @staticmethod
    def _get_saved_fig_opts_matplotlib() -> dict[str, Any]:
        return {
            # "format": "pdf",
            "bbox_inches": "tight",
            # "transparent": True,
            # "dpi": 600,
        }

    def print_steps_time(
            self,
            tag_tips: str = "Time of each step",
            tag_step_time: str = "Time of steps",
            tag_total_time: str = "Total Time",
            tag_total_iteration: str = "Total Iteration"
        ):
        """Show total time-consuming of each step
        """

        total_time = self.running_time
        total_iteration = self.kkt_iteration[-1]
        step_labels, step_time = list(self.steps_time.keys()), list(self.steps_time.values())
        sum_step_time = sum(step_time)
        sum_step_time_100_iteration = 100.0 * sum_step_time / total_iteration

        # Get the max length of step labels
        max_len = max([len(label) for label in step_labels + [tag_step_time, tag_total_time, tag_total_iteration]])

        # Cat the message
        msg_step_time = "\n".join([
            f"{step_label:<{max_len}}: {step_time:>7.2f} sec ({100.0 * step_time / total_time:5.2f}%) "
            f"({100.0 * step_time / total_iteration:<5.2f} sec/100-iterations)"
            for step_label, step_time in zip(step_labels, step_time)
        ])
        msg_total_time = \
            f"{tag_step_time.ljust(max_len)}: {sum_step_time:>7.2f} sec ({100.0 * sum_step_time / total_time:5.2f}%) ({sum_step_time_100_iteration:<5.2f} sec/100-iterations)\n"\
            f"{tag_total_time.ljust(max_len)}: {total_time:>7.2f} sec ({100.0:5.2f}%)\n"\
            f"{tag_total_iteration.ljust(max_len)}: {total_iteration:>7.0f} iterations"

        logging.info(
            f"{self._separate_symbol(tag_tips)}\n"
            f"{msg_step_time}\n"
            f"{self._sub_separate_symbol}\n"
            f"{msg_total_time}\n"
            + '=' * 60 + '\n'
        )

    def print_end_history(self):
        """Print the running history at the end
        """

        # Relative kkt errors
        kkt_errors = self.kkt_errors[-1, :]
        max_len = max([len(label) for label in self.kkt_labels])
        msg_kkt = "\n".join([f"{label:<{max_len}}: {error:>6.2e}" for error, label in zip(kkt_errors, self.kkt_labels)])

        logging.info(
            f"{self._separate_symbol("The kkt errors at end")}\n"
            f"{msg_kkt}"
        )

        # Other running history
        if self.history:
            msg_history = "\n".join([f"{key}: {history_item[-1]:.6e}" for key, history_item in self.history.items()])
            logging.info(
                f"{self._separate_symbol("Other history at end")}\n"
                f"{msg_history}"
            )




if __name__ == '__main__':
    """A simple demo to use RunningHistory class
    """

    max_iteration = 1000
    iterations = 500

    # Initialization an instance for RunningHistory
    run_history = RunHistoryManager(
        max_record_numbers=max_iteration,
        kkt_labels=["Primal Feasibility", "Dual Feasibility", "Complementary"],
        name="Alternative Direction Multiplier Method",
    )

    # Synthetic kkt errors
    rng = np.random.default_rng(42)

    noise = rng.normal(loc=0.0, scale=0.5, size=iterations)
    prim_errors = [0.4 / (n + 1) ** 2 * (1.0 + noise[n] / (n+1) ** 0.5) for n in range(iterations)]

    noise = rng.normal(loc=0.0, scale=0.5, size=iterations)
    dual_errors = [0.5 / (n + 1) ** 2 * (1.0 + noise[n] / (n+1) ** 0.5) for n in range(iterations)]

    noise = rng.normal(loc=0.0, scale=0.5, size=iterations)
    comp_errors = [0.6 / (n + 1) ** 2 * (1.0 + noise[n] / (n+1) ** 0.5) for n in range(iterations)]

    # Fake computing (to simulate ADMM)
    step1 = lambda n: time.sleep(0.001)
    step2 = lambda n: time.sleep(0.002)
    step3 = lambda n: time.sleep(0.003)

    # Collect running history
    run_history.start()
    run_history.create_tol_progress(target_tol=1e-6)
    for n in range(iterations):
        with run_history.timer(tag="Step1"):
            step1(n)

        with run_history.timer(tag="Step2"):
            step2(n)

        with run_history.timer(tag="Step3"):
            step3(n)

        kkt_errors = [prim_errors[n], dual_errors[n], comp_errors[n]]
        run_history.record(kkt_errors=kkt_errors, current_it=n)
        run_history.show_tol_progress(n, max(kkt_errors))

    # Trim running history
    run_history.end()

    # Show running history
    run_history.plot_kkt_errors(x_axis='iteration')
    run_history.plot_kkt_errors(x_axis='time')
    run_history.print_end_history()
    run_history.print_steps_time()
