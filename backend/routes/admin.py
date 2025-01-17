from fastapi import APIRouter, BackgroundTasks, Depends
from fastapi.responses import FileResponse

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
    return lowe_state_monitor


@router.get("/log")
def get_log_file():
    return FileResponse("expressinhost.log")


@router.get("/log/backup/{number}")
def get_backup_log_file(number: int):
    return FileResponse(f"expressinhost.log.{number}")
