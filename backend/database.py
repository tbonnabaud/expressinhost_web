from typing import Annotated

from fastapi import Depends, status
from fastapi.exceptions import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError, NoResultFound, SQLAlchemyError
from sqlalchemy.orm import Session, sessionmaker

from .logger import logger
from .settings import settings

DATABASE_URL = f"postgresql+psycopg://{settings.POSTGRES_USER}:{settings.POSTGRES_PASSWORD}@{settings.DB_HOST}/{settings.POSTGRES_DB}"

engine = create_engine(DATABASE_URL, echo=False, pool_size=20)
LocalSession = sessionmaker(engine)


def get_session():
    """Get a SQLAlchemy session."""
    session = LocalSession()

    try:
        yield session

    except NoResultFound as not_found_error:
        logger.error(not_found_error)
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Resource not found")

    except SQLAlchemyError as error:
        logger.error(error)
        raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, "Database error")

    finally:
        session.close()


def get_session_with_commit():
    """Get a SQLAlchemy session with commit at the end."""
    session = LocalSession()

    try:
        yield session
        session.commit()

    except IntegrityError as integrity_error:
        logger.error(integrity_error)
        session.rollback()
        raise HTTPException(status.HTTP_409_CONFLICT, "Integrity error")

    except NoResultFound as not_found_error:
        logger.error(not_found_error)
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Resource not found")

    except SQLAlchemyError as error:
        logger.error(error)
        session.rollback()
        raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, "Database error")

    finally:
        session.close()


SessionDependency = Annotated[Session, Depends(get_session)]
SessionWithCommitDependency = Annotated[Session, Depends(get_session_with_commit)]
