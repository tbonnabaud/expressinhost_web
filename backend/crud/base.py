from sqlalchemy.orm import Session


class BaseRepository:
    def __init__(self, session: Session) -> None:
        self.session = session

    def execute_with_commit(self, stmt):
        """
        This function executes a statement using a session and commits the transaction if successful,
        otherwise it rolls back the transaction.

        Args:
          stmt: The `stmt` parameter in the `execute_with_commit` method is typically a SQL statement or
        query that you want to execute using the database session. This statement can be a SELECT,
        INSERT, UPDATE, DELETE, or any other valid SQL command that you want to run against the
        database. The method

        Returns:
          The `result` variable is being returned.
        """
        try:
            result = self.session.execute(stmt)

        except:
            self.session.rollback()
            raise

        else:
            self.session.commit()
            return result
