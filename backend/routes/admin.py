from fastapi import APIRouter, Depends, BackgroundTasks

from ..authentication import check_is_admin
from ..external_db_extractors.lowe_lab import run_scraping


router = APIRouter(
    tags=["Admin"], prefix="/admin", dependencies=[Depends(check_is_admin)]
)


@router.post("/external-db/web-scraping/run")
async def run_web_scraping(background_tasks: BackgroundTasks):
    background_tasks.add_task(run_scraping)

    return {"msg": "Web scraping is running"}


@router.get("external-db/web-scraping/state")
def get_web_scraping_state():
    return {}
