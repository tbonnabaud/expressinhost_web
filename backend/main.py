from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from .database import engine
from .models import Base
from .routes import (
    codon_tables,
    codon_translations,
    results,
    tuned_sequences,
    tuning,
    users,
)

# from fastapi.middleware.cors import CORSMiddleware


Base.metadata.create_all(engine, checkfirst=True)

app = FastAPI(title="ExpressInHost")

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

api_app = FastAPI(title="ExpressInHost API")

api_app.include_router(users.router)
api_app.include_router(codon_tables.router)
api_app.include_router(codon_translations.router)
api_app.include_router(results.router)
api_app.include_router(tuned_sequences.router)
api_app.include_router(tuning.router)


@api_app.get("/")
async def api_root() -> RedirectResponse:
    return RedirectResponse(url="/api/docs")


# Serve API
app.mount("/api", api_app)

# Serve the Vue.js application
app.mount("/", StaticFiles(directory="frontend/dist", html=True), name="static")


@app.exception_handler(404)
async def root(request: Request, exc: HTTPException) -> FileResponse:
    """Catch other routes."""
    return FileResponse("frontend/dist/index.html", media_type="text/html")
