from typing import Annotated

from fastapi import Depends, status
from fastapi.exceptions import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import Session, sessionmaker

from .settings import settings

engine = create_engine(settings.DATABASE_URL, echo=False)
LocalSession = sessionmaker(engine)


def get_session():
    """Get a SQLAlchemy session."""
    session = LocalSession()

    try:
        yield session

    except SQLAlchemyError as error:
        print(error)
        raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, str(error))

    finally:
        session.close()


SessionDependency = Annotated[Session, Depends(get_session)]
