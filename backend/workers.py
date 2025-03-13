import os
import signal
from multiprocessing import Process

import redis
from rq import Worker

# Establish a connection to Redis
redis_conn = redis.Redis()
processes: list[Process] = []


def start_worker(queue_names: list[str]):
    """
    Starts an RQ worker for the specified queue names.

    Args:
        queue_names (list[str]): A list of queue names that the worker should listen to.
    """
    worker = Worker(queue_names, connection=redis_conn)
    worker.work()


def stop_workers(signal, frame):
    """
    Gracefully stops all worker processes and closes the Redis connection.

    This function is intended to be used as a signal handler for SIGINT and SIGTERM signals.

    Args:
        signal (int): The signal number that triggered this handler.
        frame (Any): The current stack frame.
    """
    print("Stopping workers...")
    for p in processes:
        p.terminate()

    for p in processes:
        p.join()

    redis_conn.close()
    print("All workers have been stopped.")


def start_and_add_process(queue_names: list[str], processes: list[Process]):
    """
    Starts a new process with the specified queue names and adds it to the list of processes.

    This function creates a new process that targets the `start_worker` function with the provided
    queue names as arguments. The process is then started and appended to the list of active processes.

    Args:
        queue_names (list[str]): A list of queue names that the new process will use.
        processes (list[Process]): A list of currently active processes to which the new process will be added.
    """
    p = Process(target=start_worker, args=(queue_names,))
    p.start()
    processes.append(p)


if __name__ == "__main__":
    # Determine the number of CPUs and set the number of workers
    num_cpus = os.cpu_count()
    num_workers = num_cpus - 1 if num_cpus > 1 else 1

    # Start a worker for web scraping and heavy processes
    # Web scraping has priority
    start_and_add_process(["web_scraping", "heavy"], processes)

    # Start workers for light processes
    for _ in range(num_workers):
        start_and_add_process(["light"], processes)

    # Set up signal handlers to stop workers gracefully on interrupt signals
    signal.signal(signal.SIGINT, stop_workers)
    signal.signal(signal.SIGTERM, stop_workers)

    # Wait for all processes to complete
    for p in processes:
        p.join()
