from fastapi import FastAPI

from .database import engine
from .models import Base

Base.metadata.create_all(engine, checkfirst=True)

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello World"}
