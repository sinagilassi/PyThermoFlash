import psutil
import os
from functools import wraps


def measure_peak_memory(func):
    """
    Decorator that prints:
    - Memory used by the function
    - System total RAM
    - System RAM usage % (before → after)
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / (1024 ** 3)  # GB

        # System memory before
        sys_mem = psutil.virtual_memory()
        sys_total_gb = sys_mem.total / (1024 ** 3)
        sys_percent_before = sys_mem.percent

        print(f"System RAM : {sys_total_gb:.2f} GB total")
        print(
            f"Before     : {mem_before:.3f} GB (process) | {sys_percent_before}% system used")

        # Run the actual function
        result = func(*args, **kwargs)

        # After measurements
        mem_after = process.memory_info().rss / (1024 ** 3)
        sys_percent_after = psutil.virtual_memory().percent

        mem_used = mem_after - mem_before

        print(
            f"After      : {mem_after:.3f} GB (process) | {sys_percent_after}% system used")
        print(f"Function used → {mem_used:.3f} GB of RAM")
        print(f"System RAM usage changed: {sys_percent_before}% → {sys_percent_after}% "
              f"({'+' if sys_percent_after > sys_percent_before else ''}"
              f"{sys_percent_after - sys_percent_before:.1f}%)")

        return result
    return wrapper


@measure_peak_memory
def create_big_list():
    a = list(range(50_000_000))   # ~400+ MB
    b = ["hello" for _ in range(10_000_000)]
    return len(a) + len(b)


# Just call it normally
# create_big_list()
