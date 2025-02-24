from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from .routes import (
    admin,
    auth,
    codon_tables,
    results,
    run_infos,
    tuned_sequences,
    tuning,
    users,
)

# from fastapi.middleware.cors import CORSMiddleware


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


api_app = FastAPI(title="ExpressInHost API", docs_url=None, redoc_url=None)
# api_app = FastAPI(title="ExpressInHost API")

api_app.include_router(auth.router)
api_app.include_router(admin.router, include_in_schema=False)
api_app.include_router(users.router)
api_app.include_router(codon_tables.router)
api_app.include_router(results.router)
api_app.include_router(tuned_sequences.router)
api_app.include_router(tuning.router)
api_app.include_router(run_infos.router)


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
