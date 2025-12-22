from fastapi import APIRouter, Depends
from fastapi.responses import PlainTextResponse, StreamingResponse

from ..authentication import check_is_admin
from ..crud.last_web_scraping import LastWebScrapingRepository
from ..database import SessionDependency
from ..external_db_extractors.lowe_lab import SOURCE, run_scraping
from ..job_manager import stream_job_state, web_scraping_queue

router = APIRouter(
    tags=["Admin"], prefix="/admin", dependencies=[Depends(check_is_admin)]
)


@router.post("/external-db/web-scraping/run")
async def run_web_scraping():
    # background_tasks.add_task(run_scraping)
    job = web_scraping_queue.enqueue(run_scraping)

    return job.id


@router.get("/external-db/web-scraping/state/{job_id}")
def stream_web_scraping_state(job_id: str):
    return StreamingResponse(stream_job_state(job_id), media_type="text/event-stream")


@router.get("/external-db/web-scraping/last-release")
def get_web_scraping_last_release(session: SessionDependency):
    return LastWebScrapingRepository(session).get_last_from(SOURCE)


@router.get("/log")
def get_log_file():
    try:
        with open("logs/expressinhost.log") as f:
            return PlainTextResponse(f.read())

    except OSError:
        return PlainTextResponse("Cannot open the log file", status_code=404)


@router.get("/log/backup/{number}")
def get_backup_log_file(number: int):
    try:
        with open(f"logs/expressinhost.log.{number}") as f:
            return PlainTextResponse(f.read())

    except OSError:
        return PlainTextResponse(
            f"Cannot open the backup log file {number}", status_code=404
        )
