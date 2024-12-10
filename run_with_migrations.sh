#!/bin/bash

# Run migrations, then add example codon tables if none exists in database, finally run the server
alembic upgrade head && python -m backend.add_example_tables && uvicorn backend.main:app --host 0.0.0.0 --workers $(nproc)