import psutil
import gc
import os

def get_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return f"Memory usage: {mem_info.rss / (1024 ** 2):.2f} MB"

    gc.collect()
