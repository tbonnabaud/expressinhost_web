import asyncio
import json

from fastapi.exceptions import HTTPException
from redis import Redis
from rq import Queue
from rq.exceptions import InvalidJobOperation, NoSuchJobError
from rq.job import Job, JobStatus

from .settings import settings

# Establish a connection to Redis
redis_conn = Redis(host=settings.REDIS_HOST, port=settings.REDIS_PORT)

# Queues
light_queue = Queue("light", connection=redis_conn)
heavy_queue = Queue("heavy", connection=redis_conn)
web_scraping_queue = Queue("web_scraping", connection=redis_conn)


def update_job_meta(job: Job, message: str, step: int, total: int | None = None):
    job.meta["message"] = message
    job.meta["step"] = step

    if total is not None:
        job.meta["total"] = total

    job.save_meta()


def get_job_state(job_id: str):
    try:
        job = Job.fetch(job_id, connection=redis_conn)
        state = {
            "status": job.get_status(),
            "step": job.meta.get("step", 0),
            "total": job.meta.get("total", 0),
        }

        if job.result:
            state["result"] = job.result

        return state

    except NoSuchJobError:
        raise HTTPException(status_code=404, detail="Job not found")


async def stream_job_state(job_id: str):
    try:
        job = Job.fetch(job_id, connection=redis_conn)

        while True:
            job.refresh()
            status = job.get_status(refresh=False)
            meta = job.get_meta(refresh=False)

            state = {
                "status": status,
                "message": meta.get("message", ""),
                "step": meta.get("step", 0),
                "total": meta.get("total", 0),
            }

            if status == JobStatus.FINISHED:
                state["message"] = "Finished."
                state["result"] = job.result
                yield json.dumps(state, default=str)
                break

            elif status in [JobStatus.FAILED, JobStatus.CANCELED]:
                state["message"] = "Failed."
                state["exc_info"] = job.exc_info
                yield json.dumps(state)
                break

            if status == JobStatus.QUEUED:
                position = job.get_position()
                state["message"] = f"Enqueued at position {position}."
                yield json.dumps(state)
                await asyncio.sleep(0.5)

            else:
                yield json.dumps(state)
                await asyncio.sleep(0.5)

    except (InvalidJobOperation, NoSuchJobError):
        yield json.dumps({"status": "not found", "message": "Not found."})
