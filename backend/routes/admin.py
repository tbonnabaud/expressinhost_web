import time

from fastapi import APIRouter, BackgroundTasks, Depends
from fastapi.responses import PlainTextResponse, StreamingResponse

from ..authentication import check_is_admin
from ..external_db_extractors.lowe_lab import lowe_state_monitor, run_scraping

router = APIRouter(
    tags=["Admin"], prefix="/admin", dependencies=[Depends(check_is_admin)]
)


@router.post("/external-db/web-scraping/run")
async def run_web_scraping(background_tasks: BackgroundTasks):
    background_tasks.add_task(run_scraping)

    return {"msg": "Web scraping is running"}


@router.get("/external-db/web-scraping/state")
def get_web_scraping_state():
    def stream():
        while True:
            yield lowe_state_monitor.model_dump_json()
            time.sleep(2)

    return StreamingResponse(stream(), media_type="text/event-stream")


@router.get("/log")
def get_log_file():
    try:
        with open("expressinhost.log") as f:
            return PlainTextResponse(f.read())

    except OSError:
        return PlainTextResponse("Cannot open the log file", status_code=404)


@router.get("/log/backup/{number}")
def get_backup_log_file(number: int):
    try:
        with open(f"expressinhost.log.{number}") as f:
            return PlainTextResponse(f.read())

    except OSError:
        return PlainTextResponse(
            f"Cannot open the backup log file {number}", status_code=404
        )
