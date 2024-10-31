# Stage 1: Build the VueJS client
FROM node:20 AS frontend-builder

# Set the working directory
WORKDIR /app/frontend

# Copy the package.json and package-lock.json files
COPY frontend/package*.json ./

# Install dependencies
RUN npm install

# Copy the rest of the frontend source code
COPY frontend/ .

# Build the frontend
RUN npm run build

# Stage 2: Final stage to serve the application
FROM python:3.12-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y dist-upgrade && \
    apt-get clean && apt-get autoremove && rm -rf /var/lib/apt/lists/*

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY requirements.txt .

RUN python -m pip install -r requirements.txt

COPY ./codon_tables ./codon_tables
COPY ./backend ./backend
COPY --from=frontend-builder /app/frontend/dist /app/frontend/dist

CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0"]