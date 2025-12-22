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
FROM python:3.13-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y dist-upgrade && apt-get install -y curl dssp && \
    apt-get clean && apt-get autoremove && rm -rf /var/lib/apt/lists/*

# Files for dssp
RUN curl -o /var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
RUN curl -o /var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
RUN curl -o /var/cache/libcifpp/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic

# Install ViennaRNA package
RUN curl -L -o /tmp/viennarna_2.7.0-1_amd64.deb "https://www.tbi.univie.ac.at/RNA/download/debian/debian_12/viennarna_2.7.0-1_amd64.deb"
RUN dpkg -i /tmp/viennarna_2.7.0-1_amd64.deb || apt-get install -y -f
RUN rm /tmp/viennarna_2.7.0-1_amd64.deb

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1
# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY requirements.txt .

RUN python -m pip install -r requirements.txt

COPY ./codon_tables ./codon_tables
COPY ./examples ./examples
COPY ./backend ./backend

# Copy migration stuff
COPY ./alembic ./alembic
COPY ./alembic.ini ./alembic.ini
COPY ./run_with_migrations.sh ./run_with_migrations.sh

COPY --from=frontend-builder /app/frontend/dist /app/frontend/dist

EXPOSE 8000

# CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0"]
CMD ["./run_with_migrations.sh"]