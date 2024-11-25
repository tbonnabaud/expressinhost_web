from typing import Annotated

from fastapi import Depends, status
from fastapi.exceptions import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm import Session, sessionmaker

from .logger import logger
from .settings import settings

DATABASE_URL = f"postgresql://{settings.POSTGRES_USER}:{settings.POSTGRES_PASSWORD}@{settings.DB_HOST}/{settings.POSTGRES_DB}"

engine = create_engine(DATABASE_URL, echo=False)
LocalSession = sessionmaker(engine)


def get_session():
    """Get a SQLAlchemy session."""
    session = LocalSession()

    try:
        yield session

    except IntegrityError as integrity_error:
        logger.error(integrity_error)
        raise HTTPException(status.HTTP_409_CONFLICT, "Integrity error")

    except SQLAlchemyError as error:
        logger.error(error)
        raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, "Database error")

    finally:
        session.close()


SessionDependency = Annotated[Session, Depends(get_session)]
