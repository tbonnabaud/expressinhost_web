from fastapi import FastAPI
from fastapi.responses import RedirectResponse

from .database import engine
from .models import Base
from .routes import codon_tables, codon_translations, results, tuned_sequences, users

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

app.include_router(codon_tables.router)
app.include_router(codon_translations.router)
app.include_router(results.router)
app.include_router(tuned_sequences.router)
app.include_router(users.router)


@app.get("/")
async def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")
