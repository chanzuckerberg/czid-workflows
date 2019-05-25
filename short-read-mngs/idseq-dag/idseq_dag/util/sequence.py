"""
  Some convenience methods to manipulate sequence objects (ex: list, range, set, string)
"""
def chunks(l, n):
    """
    Yield successive n-sized chunks from list l.
    Example:
      l = [1,2,3,4,5,6,7,8,9,10]
      c = list(chunks(l, 3))
      # c will be equal "[ [1,2,3], [4,5,6], [7,8,9], [10]]
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]
