from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .settings import settings

engine = create_engine(settings.DATABASE_URL, echo=True)
Session = sessionmaker(engine)


def get_session():
    """Get a SQLAlchemy session."""
    with Session() as session:
        yield session
