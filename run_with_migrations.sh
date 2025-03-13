#!/bin/bash

# Create logs directory
mkdir -p logs
# Run database migrations
alembic upgrade head
# Add example codon tables if none exists in database
python -m backend.add_example_tables
# Run Python RQ workers
python -m backend.workers &
# Run the periodic web scraping script in background 
python -m backend.external_db_extractors.lowe_lab &
# Run the server
uvicorn backend.main:app --host 0.0.0.0 --workers $(nproc)
