from typing import Annotated

from fastapi import Depends, Query
from sqlalchemy.orm import Session

from .database import get_session
from .schemas import FilterParams

SessionDependency = Annotated[Session, Depends(get_session)]
FilterDependency = Annotated[FilterParams, Query()]
