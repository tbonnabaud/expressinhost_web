import asyncio
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from .external_db_extractors.lowe_lab import periodic_web_scraping
from .logger import logger
from .routes import admin, auth, codon_tables, results, tuned_sequences, tuning, users

# from fastapi.middleware.cors import CORSMiddleware

background_tasks: set[asyncio.Task] = set()


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Start up
    task = asyncio.create_task(periodic_web_scraping())
    background_tasks.add(task)
    logger.info("Background periodic web scraping started.")

    yield  # Run the FastAPI application

    # Shutdown
    for task in background_tasks:
        task.cancel()

    await asyncio.gather(*background_tasks, return_exceptions=True)
    logger.info("Background tasks shut down complete.")


app = FastAPI(title="ExpressInHost", lifespan=lifespan)

# app.add_middleware(
#     CORSMiddleware,
#     allow_origins=[
#         "http://localhost",
#         "http://localhost:5173",
#         "http://localhost:8080",
#     ],
#     allow_credentials=True,
#     allow_methods=["*"],
#     allow_headers=["*"],
# )


api_app = FastAPI(title="ExpressInHost API", docs_url=None, redoc_url=None)

api_app.include_router(auth.router)
api_app.include_router(admin.router, include_in_schema=False)
api_app.include_router(users.router)
api_app.include_router(codon_tables.router)
api_app.include_router(results.router)
api_app.include_router(tuned_sequences.router)
api_app.include_router(tuning.router)


@api_app.get("/")
async def api_root() -> RedirectResponse:
    return RedirectResponse(url="/api/docs")


# Serve API
app.mount("/api", api_app)

# Serve examples files
app.mount("/examples", StaticFiles(directory="examples"), name="examples")

# Serve the Vue.js application
app.mount("/", StaticFiles(directory="frontend/dist", html=True), name="static")


@app.exception_handler(404)
async def root(request: Request, exc: HTTPException) -> FileResponse:
    """Catch other routes."""
    return FileResponse("frontend/dist/index.html", media_type="text/html")
