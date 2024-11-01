from typing import Annotated

from fastapi import Depends


class FilterParams:
    def __init__(self, offset: int | None = None, limit: int | None = None):
        self.offset = offset
        self.limit = limit


FilterParamDependency = Annotated[FilterParams, Depends(FilterParams)]
